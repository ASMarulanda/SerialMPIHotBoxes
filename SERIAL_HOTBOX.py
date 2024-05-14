
#Import required libraries 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time


nmax = 20                              # Maximum number of energy levels
N = 200                                # Number of wells
l = 10                                 # Width of the well in nm
m = 4                                  # Number of repetitions for statistics
L = np.array([2, 4, 6, 8, 10, 12])     # Array of widths in nm
Ta = np.array([10, 100, 500, 1000])    # Array of temperatures
Tf = np.array([10, 30, 50, 60, 80, 100, 200, 300, 400, 500, 600, 700, 900]) # Array of temperatures
T = np.arange(10, 1500, 50) # Array of temperatures
# L = np.arange(1, 20, 1) 
CV = np.zeros(len(Ta))   # Array to save results
EAVG = np.zeros(len(Ta))  # Array to save results


# Is it only valid in jupyter notebook? %%time

def Pot_energy(T, N, l, nmax, plot = False, live = False):
    """
    Perform a simulation to compute the potential energy of a system of N potential wells.

    Parameters:
        T (float): Temperature of the system.
        N (int): Number of potential wells in the system.
        l (float): Length of the wells.
        nmax (int): Maximum level of occupancy for a well.
        plot (bool): Whether to generate a plot of the simulation results.
        live (bool): Whether to generate an animation of the distribution of energy states.

    Returns:
        tuple: A tuple containing the average energy <E> and the heat capacity Cv.
    
    """

    # Convert well width to meters
    l = l * 1e-9

    # Constants
    kb = 8.617e-5  # Boltzmann constant in eV*K
    h = 4.11e-15   # Planck constant in eV*s
    m = 0.510      # Electron mass in MeV*c**2
    c = 3e8        # Speed of light in m/s
    k = h**2/8/l**2/m*c**2/1e6  #Energy of the ground state

    # Variables for energy and convergence
    avg_energy = 0.
    E2 = 0.
    tc = 0
    nc = 0

    # Initialize the system with random initial states
    n = [np.random.randint(1, nmax) for _ in range(N)]
    n = np.array(n)

    # Energy of the initial state
    E = [sum(n**2)*k]

    # Convergence condition
    R = True

    # Epoch initialization
    t = 1

    # Array of epochs
    y = [t]

    # Array with information about occupancy levels
    nt = []

    # Array to store cumulative energy
    EC = [sum(n**2)*k]

    start_time = time.time()  # Record start time

    # Main simulation loop
    while t < 1500:
        t += 1  # Increase epoch

        Ee = 0  # Initialize energy in epoch t

        if live and t<100:
          # Generate live animation if enabled
            plt.clf()
            plt.hist(n, bins=15, range=(1,15), label='t = %s'%(t))
            plt.title('Histogram of Energy Levels Occupation at T = %s' %(T))
            plt.xlabel('Energy Level')
            plt.ylabel('Frequency')
            plt.ylim(0,50)
            plt.legend()
            plt.pause(0.5)
            #plt.savefig(f'histogram_animation_frame_{t}.png')

        # Perform random walk for each potential well
        for i in range(N):
            r = np.random.uniform()

            # Randomly choose to decrease or increase occupancy
            if r > 0.5 and n[i] > 1:
                nc = n[i] - 1
            else:
                nc = n[i] + 1

            # Calculate energy changes
            Ec = k*nc**2
            Ei = k*n[i]**2
            dE = Ec - Ei

            # Accept or reject the transition based on temperature
            if dE <= 0:
                n[i] = nc
            else:
                r = np.random.uniform()
                W = np.exp(-dE/T/kb)
                if r <= W:
                    n[i] = nc

            # Update energy in epoch t
            Ei = k*n[i]**2
            Ee += Ei

        # Calculate cumulative energy
        Ecum = sum(E)/t
        EC.append(Ecum)

        # Check for convergence
        if t > 100:
            nt.append(n[i])

        if t == 1499:
          e2 = np.array(E[-1000:])
          E2 = (sum(e2**2))/(len(e2))
          avg_energy = np.mean(e2)

        # Store energy for epoch t
        E.append(Ee)

        # Calculate cumulative energy
        Ecum = sum(E)/t
        EC.append(Ecum)

        # Append epoch to array
        y.append(t)


    elapsed_time = time.time() - start_time  # Calculate elapsed time
    elapsed_times.append(elapsed_time)       # Store elapsed time

    # Cv calculation
    cv = (E2-avg_energy**2)/kb/T**2

    # Plot results if required. 
    if plot:
        fig = plt.figure(figsize=(10,4))
        fig.suptitle('Evolución de la energía de un sistema de 20 pozos infinitos de potencial \ny la distribución de los niveles de energía durante esta')

        # Plot energy evolution
        plt.subplot(1, 2, 1)
        plt.plot(y, E, label=f'Temperature: {T}')
        plt.xlabel('Epoch')
        plt.ylabel('Total energy (eV)')
        plt.title('Total energy evolution')
        plt.legend(loc='lower left')
        #plt.text(0.07, 0.9, convergence_info, transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.5))

        # Plot histogram of nt values
        plt.subplot(1, 2, 2)
        nt = np.array(nt)
        nt = nt - 1 #Redifine the ground level from 1 to 0 in order to make the adjustment
        counts, bins, _ = plt.hist(nt, bins=20, range=(0,20), density=True)
        bins = bins[:-1]

        plt.title('Histogram of Energy Levels')
        plt.xlabel('Energy Level')
        plt.ylabel('Frequency')
        plt.tight_layout()

        def boltz(E,a,b):
          return np.exp(-E*a/kb/T)*E*k/b

        #Make a curve fitting of the histogram to a Maxwell-Boltzman distribution
        resul, err = curve_fit(boltz, bins, counts, p0 = (k, 5 ) )
        b = np.arange(0,20,0.1)
        r=sum((counts-boltz(bins,*resul))**2)/(len(counts)-2)
        v=np.std(counts)**2
        R=1-(r/v)
        plt.plot(b, boltz(b, *resul ))
        plt.text(0.7,0.8, "R^2 = %.2f" %(R),  transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.5))


        plt.show()

    return avg_energy, cv

def Potential_reservoir(Ta, N, l, nmax, plot=False, live=False):
    """
    Perform a simulation to compute the total energy of an array of N potential
    wells when the system reach the equilibrium with a thermal reservoir at a temperature T.

    Parameters:
        Ta (float, list, np.ndarray): Temperature or array of temperatures.
        N (int): Number of potential wells in the system.
        l (float, list, np.ndarray): Length of the wells or array of lengths in nm.
        nmax (int): Maximum level of occupancy for a well.
        plot (bool): Whether to generate a plot of the simulation results.
        live (bool): Whether to generate an animation of the distribution of energy states.

    Returns:
        tuple: A tuple containing arrays of average energy and heat capacity.
    """

    # Lists to store average energy and heat capacity
    E_avg = []
    Cv = []

    # Check if Ta is an array of temperatures
    if isinstance(Ta, (list, np.ndarray)):
        for i in Ta:
            E, C = Pot_energy(i, N, l, nmax, plot, live)
            E_avg.append(E)
            Cv.append(C) #Array of the heat capacity values: the function does not return this array beacuse the result is not convincing
        return E_avg, Cv

    # Check if l is an array of lengths
    elif isinstance(l, (list, np.ndarray)):
        for i in l:
            E, C = Pot_energy(Ta, N, i, nmax, plot, live)
            E_avg.append(E)
            Cv.append(C)  #Array of the heat capacity values: the function does not return this array beacuse the result is not convincing
        return E_avg, Cv

    # Otherwise, perform simulation with single temperature and length
    else:
        E_avg, Cv = Pot_energy(Ta, N, l, nmax, plot, live)
        return E_avg , Cv


#Definition of the functions to fit
def line(x,a,b):
  return a*x + b

def e(x,a,b,c):
  return np.exp(-x*a)*b+c


elapsed_times = []  # List to store elapsed times

# Lists to store average energy and heat capacity
CV = []
EAVG = []

start_time = time.time()  # Record start time

# For the specified number of repetitions (m), each iteration calculates the average energy (E_avg) and heat capacity (Cv) 
# for a system with N wells, each of width l, at an initial temperature Tf in the thermal reservoir. This process aims to improve 
# the final results by reducing standard deviation and enhancing efficiency through multiple runs and averaging of results.
for i in range(m):       
    E_avg, Cv = Potential_reservoir(Tf, N, l, nmax)
    CV.append(Cv)
    EAVG.append(E_avg)

# For serial computing, graphical processes can be invoked.
fig = plt.figure(figsize=(10,4))

# Loop over the results
for i in range(len(EAVG)):
    # Plot total energy vs temperature
    plt.subplot(1, 2, 1)
    plt.plot(Tf, EAVG[i], label=r'run: {}'.format(i+1), color=['darkblue','limegreen','orange', 'darkred'][i])
    # Set title and labels
    plt.title('Total energy vs Temperature')
    plt.ylabel('Total average energy (eV)')
    plt.xlabel('Temperature')
    plt.legend()

    # Plot heat capacity vs temperature
    plt.subplot(1, 2, 2)
    plt.plot(Tf, CV[i], label=r'run: {}'.format(i+1), color=['darkblue','limegreen','orange', 'darkred'][i])
    # Set title and labels
    plt.xlabel('Temperature')
    plt.ylabel('Cv (eV/K)')
    plt.title('Heat capacity vs temperature')

# Display legends and adjust layout
plt.legend()
plt.tight_layout()
plt.show()

# Calculate and display total elapsed time
elapsed_time = time.time() - start_time
print(f"Total time: {elapsed_time}")