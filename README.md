# Monte Carlo Simulation of Particle Distribution in Potential Wells with Thermal Reservoir

We present the adaptation of a Monte Carlo-based simulation for obtaining the energy distribution of particles within an array of potential wells in contact with a thermal reservoir. 
The main objective of this project is to investigate the efficiency of the simulation when executed using both single-processor serial and parallel processing via SSH cluster connection. By dividing the algorithm into parallelized tasks, we explore the benefits of parallel computing, focusing on computational time and resource utilization.  The final results use simultaneous runs to display average energy behavior according to temperature increase, the same as the system's heat capacity.

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/6417c58f-6d63-45cb-9f46-21eb815a68d7" width="400">


Traditionally, programming favored serial computation, where problems are tackled sequentially on a single processor. However, with growing computing demands, parallel computing has become crucial. It utilizes multiple processors to solve problems concurrently, dividing tasks into manageable parts executed simultaneously. This approach enhances performance and reduces execution time for large-scale computations. Effective coordination mechanisms ensure seamless collaboration among processors, optimizing overall efficiency.

##  Monte Carlo-Metropolis Simulation 

Monte Carlo simulations employ randomness and probabilistic sampling to tackle problems where deterministic solutions are challenging to obtain. A notable example is the Metropolis algorithm, widely utilized in statistical physics. This algorithm generates energy values consistent with a Maxwell-Boltzmann distribution to derive expected values for observables. It employs Markov Chain Monte Carlo (MCMC) at its core, exploring a system's phase space and sampling configurations based on their Boltzmann weights. By sequentially proposing changes and accepting or rejecting them probabilistically, MCMC ensures equilibrium is reached and samples configurations according to their equilibrium probabilities, adhering to detailed balance conditions. The likelihood of being in a particular state is determined by the system's energy and temperature.

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
        Tuple: Array or integer representing the average energy. The format of the tuple varies depending on the type of the Ta parameter: if Ta is an integer, both average energy and heat capacity are single values; if Ta is a list, each element of the tuple corresponds to an array of values, and similarly for the l parameter.
    

This is the main function for the user, from it the user can obtain how is the behaviour of the system energy respect to the different parameters.


# Heat Capacity 
The heat capacity of a system can be computed from its energy behavior when it is in thermal equilibrium, where $k_b$ is the Boltzmann constant, as:

$C_v=\frac{\left\langle E^2\right\rangle - \langle E\rangle^2}{k_b T^2}$.

# Results - Graphics
![00SERIALRESULTS](https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/7c87ce1e-0009-42b5-b408-3cda69b97bb2)


