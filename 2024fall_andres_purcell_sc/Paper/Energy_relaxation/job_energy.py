import two_sites_H2 as RETH2
import numpy as np
import sys

size = (10, 10) # TLS polarization box
c1_l, c2_l = np.sqrt(1.0), np.sqrt(0.0) #TLS initial state

tend = 40000 # Total simulation time
pml_length = 2.0
mu12 = 0.03 # Dipole transition moment
#omega_0 = 0.10 # TLS transition energy frequency

print_every = 1 # Save data each time step

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
    omega_0 = 0.10,
    counter=0,
    gauss_freq=1.0,
    outputName="tmp.txt",
    outputName2="tmp.txt"
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
        gauss_freq=gauss_freq,
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
    sys.output_EM_energy(outputName2)

if __name__ == "__main__":

    sigmas = np.array([0,5,10,15,20,25,30,35])*0.001

    # Up and LP frequencies for Lorentzian_Gamma = 1e-5
    UP_Freq = [0.1, 0.10260863, 0.10370163, 0.10457063, 0.10532663, 0.10597163, 0.10658863, 0.10714863]
    LP_Freq = [0.1, 0.09751938, 0.09656638, 0.09586538, 0.09527638, 0.09477238, 0.09429538, 0.09390338]

    # Up and LP frequencies for Lorentzian_Gamma = 1e-4
    UP_Freq2 = [0.1000,0.10256281407035175,0.10369346733668341,0.10457286432160805,0.10532663316582914,0.10595477386934674,0.10658291457286433,0.10721105527638192]
    LP_Freq2 = [0.1000,0.09766331658291458,0.09665829145728644,0.09590452261306533,0.09540201005025126,0.09489949748743719,0.09439698492462312,0.09402010050251255]

    # Up and LP frequencies for Lorentzian_Gamma = 1e-3
    UP_Freq3 = [0.1, 0.10268844, 0.10383567, 0.10471667, 0.10582915, 0.1060804 , 0.10670854, 0.10733668]
    LP_Freq3 = [0.1, 0.09778894, 0.09681868, 0.09609768, 0.09653266, 0.0948995 , 0.09452261, 0.09414573]

    counter = int(sys.argv[-1])
    print("counter from bash is", counter)

    R = 6
    n1 = 2.0  # mirrors
    n2 = 1.0  # air
    rescaling = 10.0  # 2 times the distance between the first mirrors of the cavity
    epsilon = 1.0  # Lorentzian medium epsilon
    Medium_Frequency = 0.1 / (2 * np.pi)
    Lorentzian_Gamma = 1e-5
    Lorentzian_Sigma = 0.02

    #for n_sigmas in sigmas:
    n_sigmas = len(sigmas)
    print("the total size of the sigmas array is", n_sigmas)
    if counter < n_sigmas:
        sigma_value = sigmas[counter]
        Freq_value = UP_Freq[counter]
        print("the sigma value now is", sigma_value)
        perform_Eh_RETH2(
            R,
            rescaling,
            n1,
            n2,
            epsilon,
            pml_length,
            Medium_Frequency,
            Lorentzian_Gamma,
            sigma_value,
            0.10532663316582914, #TLS frequency
            counter,
            Freq_value,
            outputName="data/EM_Traj2_UP_Sigma_%3f.txt" %(sigma_value),
            outputName2="data/EM_Energy2_UP_Sigma_%3f.txt" %(sigma_value)
        )

# suppression frequency at 0.10532663316582914 & LP frequency at 0.09527638190954774 for Lorentzian_Gamma = 1e-5
# suppression frequency at 0.10532663316582914 & LP frequency at 0.09549291995925126 for Lorentzian_Gamma = 1e-4
# suppression frequency at 0.10532663316582914 & LP frequency at 0.09653266331658292 for Lorentzian_Gamma = 1e-3

