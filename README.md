# Paralell Adapting of Monte Carlo Simulation for Particle Distribution in Potential Wells with Thermal Reservoir

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


### Serial Results - Graphics

This initial version of the code operates sequentially, simulating the behavior of the system at each temperature one after the other. Here, the code is run four times. 

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/5f506694-4812-4f0b-b720-7f3e0e36ec21" width="400">


This graph depicts the total average energy of the system at equilibrium with the temperature of the thermal bath for four separate simulation runs.

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/dea1e137-386a-4794-a70b-bbed11ecbd98" width="400">

This visualization illustrates the heat capacity of the system at equilibrium for the temperature of the thermal bath for four distinct simulation runs. 

The magic command %%time provides essential metrics for evaluating code performance. User Time reflects the total CPU time dedicated to executing user instructions, encompassing the code and function calls. System Time quantifies CPU time allocated to executing system-related tasks, such as handling inputs/outputs. Total Time aggregates User Time and System Time, offering a comprehensive measure of overall CPU time consumption. Wall Time represents the actual elapsed time from the start to the end of script execution, including waiting periods for resources and libraries. 

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/4c893418-8fa7-403c-a790-905994ce9f32" width="400">


#  Paralell Approach 

Optimizing the performance of computing units is fundamental in parallel computing. Eliminating graphical processes, known for consuming substantial RAM, can enhance performance significantly. However, this adjustment may lead to differences in simulation output compared to previous versions. 
The code utilizes MPI's distributed computing approach to distribute simulation runs among processes, each managing simulations for assigned temperatures independently. Results are stored locally and later compiled into a (.csv) file for analysis. This strategy primarily aims to reduce computing time compared to serial computation.

On a cluster with 8 CPUs, each managing 4 simulation runs, the total computing time was notably reduced. Remarkably, results from parallel computation closely matched those from serial computation, accompanied by error bars.

### Paralell Results - Graphics
<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/6123cd9e-d820-436f-8d5f-8e6b438df760" width="400">

The linear progression of total energy with increasing temperature aligns with theoretical expectations, indicating the system's sensitivity to temperature variations. Additionally, the smooth curve observed in the heat capacity graph presents an averaged fit capturing the response of system particles to energy transitions between levels. This representation offers deeper insight into the system's dynamic behavior under varying temperature conditions, confirming the accuracy and reliability of the parallel simulation.

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/6ef9d168-7cca-464c-aa28-9104ff5428de" width="400">


### Error Analysis 
In the results, the error must be accounted for. Each data point on the graphs is accompanied by an associated error bar, representing the standard deviation calculated from the resulting values obtained on each processor for the multiple simulations.

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/ef447af3-bf54-4b95-9d07-41d0eb982e49" width="400">

The observed error in the heat capacity shows a noticeable upward trend with increasing temperature. This trend can be primarily attributed to the stochastic nature of the computed values, leading to greater uncertainty in the heat capacity estimation, especially at higher temperatures.

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/e0632778-fbb1-469d-b6fc-220e99a64ce7" width="400">

On this result, the error demonstrates a linear increase with temperature. Despite this, the error propagation remains relatively minimal compared to the uncertainty in heat capacity, particularly at high temperatures. This indicates that the total energy measurements provide highly reliable outcomes, with the averaging across multiple calculations facilitating precise forecasts of the system's dynamics.


