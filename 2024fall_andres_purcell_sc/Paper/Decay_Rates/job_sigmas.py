import two_sites_H2 as RETH2
import numpy as np
import sys

size = (10, 10) # TLS polarization box
c1_l, c2_l = np.sqrt(0.9), np.sqrt(0.1) # TLS initial state

tend = 100000 # Total simulation time
pml_length = 2.0
mu12 = 0.03 # Dipole transition moment
omega_0 = 0.09653266331658292 # TLS transition energy frequency

print_every = 30 # Save data each 30th time step

def perform_Eh_RETH2(
    R=6,
    rescaling=10.0,
    n1=2.0,
    n2=1.0,
    epsilon=1,
    pml_length=2.0,
    Medium_Frequency=1.0,
    Lorentzian_Gamma=1e-5,
    Lorentzian_Sigma=0.005,
    counter=0,
    outputName="tmp.txt"
):
    # Firstly, we perform a single Ehrenfest dynamics
    sys = RETH2.TwoSites(
        separation=R,
        omega_0=omega_0,
        mu12=mu12,
        c1_l=c1_l,
        c2_l=c2_l,
        rescaling=rescaling,
        size=size,
        n1=n1,
        n2=n2,
        pml_length=pml_length,
        print_every=print_every,
    )
    # Secondly, we propagate the electrodynamics of the system and/or integrate the total EM field
    sys.integrate(
        tmax=tend,
        output_field=False,
        epsilon=epsilon,
        Medium_Frequency=Medium_Frequency,
        Lorentzian_Gamma=Lorentzian_Gamma,
        Lorentzian_Sigma=Lorentzian_Sigma,
    )
    sys.output_traj(outputName)

if __name__ == "__main__":
    sigmas = np.array([0,1,5,10,30,50,100,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500])*1e-5

    counter = int(sys.argv[-1])
    print("counter from bash is", counter)

    R = 6
    n1 = 2.0  # mirrors
    n2 = 1.0  # air
    rescaling = 10.0  # 2 times the distance between the first mirrors of the cavity
    epsilon = 1.0  # Lorentzian medium epsilon
    Medium_Frequency = 0.1 / (2 * np.pi)
    Lorentzian_Gamma = 1e-3
    Lorentzian_Sigma = 0.02

    #for n_sigmas in sigmas:
    n_sigmas = len(sigmas)
    print("the total size of the sigmas array is", n_sigmas)
    if counter < n_sigmas:
        LSigma = sigmas[counter]
        print("the sigma value now is", LSigma)
        perform_Eh_RETH2(
            R,
            rescaling,
            n1,
            n2,
            epsilon,
            pml_length,
            Medium_Frequency,
            Lorentzian_Gamma,
            LSigma,
            counter,
            outputName="data/trajEh_LP_DECAY_C&LM_GAMMAe-3_SIGMA_%5f.txt" %LSigma
        )

