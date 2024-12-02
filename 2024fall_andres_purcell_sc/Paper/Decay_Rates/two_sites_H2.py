"""
Here, the two sites are in diabatic basis, and there is no direct dipole-dipole
interactions and we propagate the total E field.

However, we may add the +R corrections, which are purely transverse fields for
the EM field and the common +R Dissipation for the sites

The implementation for this can be very simple
"""

import meep as mp
import numpy as np
import matplotlib.pyplot as plt
from two_level_system import TLS
from math import exp


class TwoSites:
    def __init__(
        self,
        omega_0=0.25,
        mu12=0.25,
        sigma=1.0,
        c1_l=np.sqrt(0.5),
        c2_l=np.sqrt(0.5),
        c1_r=np.sqrt(1),
        c2_r=np.sqrt(0),
        rescaling=10.0,
        pml_length=2.0,
        n1=2.0,
        n2=1.0,
        size=(10.0, 10.0),
        resolution=40,
        print_every=30,
        phi_l=0.0,
        phi_r=0.0,
        gauss_freq = 1.0,
    ):

        # Initialize TLS
        self.L = TLS(
            omega_0=omega_0,
            mu12=mu12,
            sigma=sigma,
            c1_0=c1_l,
            c2_0=c2_l,
            center=(0, 0),
            size=size,
            dx=self.dx,
            dt=self.dt,
        )

        # parameters for simulation
        self.courant = 0.5
        self.rescaling = rescaling * 2.0 * np.pi
        self.resolution = resolution / self.rescaling
        self.dx = 1.0 / self.resolution
        self.dt = self.courant * self.dx
        self.print_every = print_every
        self.t = 0.0
        self.energy = 0.0
        self.count = 0
        self.new_x_range = 0
        self.frequency = (gauss_freq) / self.rescaling 
        self.frequencyWidth = 0.05 / self.rescaling
        self.size = size
        self.numberFrequencies = 200

        # prepare the sim variable
        self.sim = []
        self.geometry = []
        self.cell = mp.Vector3()
        self.pml_layers = []
        self.vals = []

        # save data
        self.t_lst = []
        self.p2L_lst = []

        # prepare the geometry
        self.pmlThickness = pml_length*self.rescaling
        self.sources = []
        self.incidentFluxToSubtract = []
        self.t1 = 0.125 * self.rescaling
        self.t2 = 0.25 * self.rescaling
        self.n1 = n1
        self.n2 = n2

        self.layerIndexes = np.array(
            [
                1.0,
                self.n1,
                self.n2,
                self.n1,
                self.n2,
                self.n1,
                self.n2,
                self.n1,
                self.n2,
                self.n1,
                1,
                self.n1,
                self.n2,
                self.n1,
                self.n2,
                self.n1,
                self.n2,
                self.n1,
                self.n2,
                self.n1,
                1.0,
            ]
        )
        self.layerThicknesses = np.array(
            [
                5*self.rescaling,
                self.t1,
                self.t2,
                self.t1,
                self.t2,
                self.t1,
                self.t2,
                self.t1,
                self.t2,
                self.t1,
                0.5 * self.rescaling,
                self.t1,
                self.t2,
                self.t1,
                self.t2,
                self.t1,
                self.t2,
                self.t1,
                self.t2,
                self.t1,
                5*self.rescaling,
            ]
        )
        self.layerThicknesses[0] += self.pmlThickness 
        self.layerThicknesses[-1] += self.pmlThickness 
        self.length = np.sum(self.layerThicknesses)
        self.layerCenters = np.cumsum(self.layerThicknesses) - self.layerThicknesses / 2
        self.layerCenters = self.layerCenters - self.length / 2
        self.sourceLocation = mp.Vector3(
            self.layerCenters[0] - self.layerThicknesses[0] / 4, 0, 0
        )
        self.transmissionMonitorLocation = mp.Vector3(
            self.layerCenters[-1] - self.pmlThickness / 2, 0, 0
        )
        self.reflectionMonitorLocation = mp.Vector3(
            self.layerCenters[0] + self.layerThicknesses[0] / 4, 0, 0
        )
        self.pmlLayers = [mp.PML(thickness=self.pmlThickness)]
        self.incidentRegion = []
        self.length
        self.cellSize = mp.Vector3(
            self.length, 3/5*self.length, 0
        )

        # Gaussian Source
        self.sources_gaussian = [
            mp.Source(
                mp.GaussianSource(frequency=self.frequency, fwidth=self.frequencyWidth),
                component=mp.Ez,
                center=self.sourceLocation,
                size=mp.Vector3(
                    0,
                    3/5*self.length-2*self.pmlThickness,
                    self.length,
                )
            )
        ]

    def propogate_Eh(self, int_EP_L):
        # For Hamiltonian II, everything is very simple
        self.L.propogate_Eh(int_ep=int_EP_L)
        self.t += self.dt
        self.count += 1

    def _step_function(self):
        def __step_func__(sim):
            if self.count % self.print_every == 0:
                self.t_lst.append(self.t)
                self.p2L_lst.append(self.L.rho[1, 1])
            # Firstly let's calculate the DP integration by MPI interface
            int_EP_L = sim.integrate_field_function(
                [mp.Ez],
                lambda R, ez: self.L._prefactor_polarization
                * exp(
                    -((R.x - self.L.center.x) ** 2 + (R.y - self.L.center.y) ** 2)
                    / (2.0 * self.L.sigma**2)
                )
                * ez,
                mp.Volume(size=self.L.size, center=self.L.center),
            )
            # Secondly, let's propogate the rho
            self.propogate_Eh(int_EP_L)
            # Thirdly, update the sources
            sim.change_sources(self.L.sources)

        return __step_func__

    def geo(
        self,
        epsilon=1.0,
        Medium_Frequency=1.0,
        Lorentzian_Gamma=1e-5,
        Lorentzian_Sigma=0.005,
    ):
        # Let's create the mirror geometry
        susceptibilities = [
            mp.LorentzianSusceptibility(
                frequency=Medium_Frequency,
                gamma=Lorentzian_Gamma,
                sigma=Lorentzian_Sigma,
            )
        ]
        Lorentzian_material = mp.Medium(
            epsilon=epsilon, E_susceptibilities=susceptibilities
        )
        idx_middle = 10
        self.geometry = [
            mp.Block(
                mp.Vector3(self.layerThicknesses[i], mp.inf, mp.inf),
                center=mp.Vector3(self.layerCenters[i], 0, 0),
                material=mp.Medium(index=self.layerIndexes[i]),
            )
            for i in range(self.layerThicknesses.size)
        ]
        self.geometry[idx_middle] = mp.Block(
            mp.Vector3(self.layerThicknesses[idx_middle], mp.inf, mp.inf),
            center=mp.Vector3(self.layerCenters[idx_middle], 0, 0),
            material=Lorentzian_material,
        )

    def integrate(
        self,
        tmax=600,
        nframes=100,
        epsilon=1.0,
        Medium_Frequency=1.0,
        Lorentzian_Gamma=1e-5,
        Lorentzian_Sigma=0.005,
        output_field=True,
    ):
        self.geo(
            epsilon,
            Medium_Frequency,
            Lorentzian_Gamma,
            Lorentzian_Sigma,
        )
        self.sim = mp.Simulation(
            cell_size=self.cellSize,
            boundary_layers=self.pmlLayers,
            geometry=self.geometry,
            sources=self.L.sources,
            resolution=self.resolution,
            Courant=self.dt / self.dx,
        )
        self.sim.run(self._step_function(), until=tmax)

    def output_traj(self, traj_filename="traj.txt"):
        # output data finnally
        with open(traj_filename, "w") as handle:
            for t, rho1 in zip(self.t_lst, self.p2L_lst):
                handle.write("%.4E %.8E \n" % (t, np.abs(rho1)))
