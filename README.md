# Monte Carlo Simulation of Particle Distribution in Potential Wells with Thermal Reservoir

We introduce a Monte Carlo-based simulation for analyzing the energy distribution of particles within an array of potential wells coupled with a thermal reservoir. The primary objective of this project is to assess the simulation's efficiency through both single-processor serial and parallel processing via SSH cluster connection. By parallelizing the algorithm, we aim to explore the advantages of parallel computing, particularly focusing on computational time and resource utilization. The final results encompass simultaneous runs illustrating the average energy behavior with increasing temperature, resembling the system's heat capacity results as well. 

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/6417c58f-6d63-45cb-9f46-21eb815a68d7" width="400">

In traditional programming, serial computation sequentially tackles problems on a single processor. However, with escalating computing demands, parallel computing has become indispensable. Utilizing multiple processors concurrently divides tasks into manageable parts, enhancing performance and reducing execution time for large-scale computations. Effective coordination mechanisms ensure seamless collaboration among processors, optimizing overall efficiency.

###  Monte Carlo-Metropolis Simulation 

Monte Carlo simulations employ randomness and probabilistic sampling to address problems with challenging deterministic solutions. The Metropolis algorithm, a prominent example, generates energy values consistent with a Maxwell-Boltzmann distribution, deriving expected values for observables. Employing Markov Chain Monte Carlo (MCMC), the algorithm explores a system's phase space, sampling configurations based on their Boltzmann weights. MCMC ensures equilibrium is reached, adhering to detailed balance conditions, and sample configurations according to their equilibrium probabilities, dependent on the system's energy and temperature.


#  Serial Approach 

The code has the function:

Potential_reservoir(Ta, N, l, nmax, plot=False, live=False)

    Perform a simulation to compute the total energy of an array of N potential wells when the system reach the equilibrium with a thermal reservoir at a temperature T. The functions make use of the function Pot_energy(T, N, l, nmax, plot, live).

    Parameters:
        Ta (float, list, np.ndarray): Temperature or array of temperatures.
        N (int): Number of potential wells in the system.
        l (float, list, np.ndarray): Length of the wells or array of lengths in nm.
        nmax (int): Maximum level of occupancy for a well.
        plot (bool): Whether to generate a plot of the simulation results. Total energy vs epoch and a histogram of the occupation distribution of the energy states. Default: False
        live (bool): Whether to generate an animation of the distribution of energy states. Default: False

    Returns:
        Tuple: Array or integer representing the average energy. The format of the tuple varies depending on          the type of the Ta parameter: if Ta is an integer, both average energy and heat capacity are single           values; if Ta is a list, each element of the tuple corresponds to an array of values, and                     similarly for the l parameter.
    
This function simulates the total energy of an array of potential wells in thermal equilibrium with a reservoir at a specified temperature. The user can input parameters such as temperature (Ta), number of wells (N), well size (l), and maximum occupancy level (nmax). The function generates an array representing the average energy. Plotting and animation options (plot and live) are available but disabled in this context. The function internally calls Pot_energy to specify system parameters and compute the potential energy.

The heat capacity (Cv) of the system can be calculated from its energy behavior in thermal equilibrium, where $k_b$ is the Boltzmann constant, as:

$C_v=\frac{\left\langle E^2\right\rangle - \langle E\rangle^2}{k_b T^2}$.


### Seriall Approach Results - Graphics

![00SERIALRESULTS](https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/7c87ce1e-0009-42b5-b408-3cda69b97bb2)


