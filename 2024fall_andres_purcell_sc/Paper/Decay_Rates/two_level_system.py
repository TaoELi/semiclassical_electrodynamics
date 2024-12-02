'''
Here we will construct the python class for the TLS.
Further, we will construct class to propograte the EM field + TLS(s) together.
'''

import meep as mp
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
from math import exp
from transverse_vector_field import calc_transverse_components

class TLS:
    def __init__(
        self,
        omega_0 = 0.10,
        mu12 = 0.10,
        sigma = 1.0,
        c1_0 = np.sqrt(0.50),
        c2_0 = np.sqrt(0.50),
        center = (.0, .0),
        size = (4., 4.),
        dx = 0.4,
        dt = 0.2,
        print_every=30,
        Phi = 0.0,
        ):
        self.omega_0, self.mu12, self.sigma = omega_0, mu12, sigma
        self.kFGR = self.omega_0**3.0 * self.mu12**2.0 / 3.0 / np.pi
        self.center = mp.Vector3(*center)
        self.size = mp.Vector3(*size)
        self._prefactor_polarization = 1.0 / (2.0 * np.pi)**(1.5) / self.sigma**2 * self.mu12 * (2*np.pi)**(0.5)
        # Parameters for simulation
        self.dx = dx
        self.dt = dt
        self.t = 0.0
        self.count = 0
        self.print_every = print_every

        # now let's construct the sources
        pz = calc_transverse_components(size=size, dx=self.dx, sigma=self.sigma, mu12=self.mu12)
        self.sources = [mp.Source(mp.CustomSource(src_func=lambda t: 1.0),
                            component=mp.Ez,
                            center=self.center,
                            size=self.size,
                            amplitude=0.0,
                            amp_data=pz
                            ) ]
        
        # Parameters for Hamiltonian
        self.Hs = np.matrix([[0, 0],[0, self.omega_0]], dtype=np.complex128)
        self.SIGMAX = np.matrix([[0,1],[1,0]], dtype=np.complex128)
        self.H = self.Hs
        self.expHs = expm(-1j * self.dt * self.Hs / 2.0)
        # Initilize the electronic population
        self.C = np.matrix([[c1_0],[c2_0]], dtype=np.complex128)
        self.rho = np.dot(self.C, self.C.conj().transpose() )
        self.U = np.matrix([[1.0, 0],[0.0, 1.0]], dtype=np.complex128)
        # prepare the sim variable
        self.sim = []
        self.cell = mp.Vector3()
        self.pml_layers = []
        self.vals = []
        self.t_lst = []
        self.p2_lst = []
        self.rho12_lst = []
        # +R parameters
        self.phi = Phi
        self.amp_R = 0.0

    def propogate_Eh(self, int_ep):
        #''' Construct the time propagator U'''
        self.U = np.dot(self.expHs, np.dot(expm( (1j * self.dt * int_ep) * self.SIGMAX), self.expHs))
        self.rho = np.dot(np.dot(self.U, self.rho), self.U.conj().transpose() )
        #''' Calculate the Ehrenfest amplitude '''
        amp = -2.0 * self.omega_0 * np.imag(self.rho[0,1])
        for s in self.sources:
            s.amplitude = amp
        # Lastly lets update time and count
        self.t += self.dt
        self.count += 1


    def _step_function(self):
        def __step_func__(sim):
            # store to data
            if self.count % self.print_every is 0:
                self.t_lst.append(self.t)
                self.p2_lst.append(np.abs(self.rho[1,1]))
            # Firstly let's calculate the EP integration by MPI interface
            int_ep = sim.integrate_field_function([mp.Ez],
                        lambda R, ez : self._prefactor_polarization * exp(-(R.x**2 + R.y**2) / (2.0 * self.sigma**2)) * (ez),
                        mp.Volume(size=self.size, center=self.center))
            # Second propogate electronic degrees of freedom using the existing coupling
            self.propogate_Eh(int_ep = int_ep)
            # thirdly we update the sources
            sim.change_sources(self.sources)
        return __step_func__

    def output_traj(self, traj_filename="traj.txt"):
    # output data finnally
        with open(traj_filename, 'w') as handle:
            for t, rho in zip(self.t_lst, self.p2_lst):
                handle.write("%.4E %.8E\n" %(t, np.abs(rho)) )
