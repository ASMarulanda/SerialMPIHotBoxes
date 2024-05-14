# Monte Carlo Simulation of Particle Distribution in Potential Wells with Thermal Reservoir

We present the adaptation of a Monte Carlo-based simulation for obtaining the energy distribution of particles within an array of potential wells in contact with a thermal reservoir. 
The main objective of this project is to investigate the efficiency of the simulation when executed using both single-processor serial and parallel processing via SSH cluster connection. By dividing the algorithm into parallelized tasks, we explore the benefits of parallel computing, focusing on computational time and resource utilization.  The final results use simultaneous runs to display average energy behavior according to temperature increase, the same as the system's heat capacity.

<img src="https://github.com/ASMarulanda/SerialMPIHotBoxes/assets/123122569/6417c58f-6d63-45cb-9f46-21eb815a68d7" width="400">


Traditionally, programming favored serial computation, where problems are tackled sequentially on a single processor. However, with growing computing demands, parallel computing has become crucial. It utilizes multiple processors to solve problems concurrently, dividing tasks into manageable parts executed simultaneously. This approach enhances performance and reduces execution time for large-scale computations. Effective coordination mechanisms ensure seamless collaboration among processors, optimizing overall efficiency.

##  Monte Carlo-Metropolis Simulation 

Monte Carlo simulations employ randomness and probabilistic sampling to tackle problems where deterministic solutions are challenging to obtain. A notable example is the Metropolis algorithm, widely utilized in statistical physics. This algorithm generates energy values consistent with a Maxwell-Boltzmann distribution to derive expected values for observables. It employs Markov Chain Monte Carlo (MCMC) at its core, exploring a system's phase space and sampling configurations based on their Boltzmann weights. By sequentially proposing changes and accepting or rejecting them probabilistically, MCMC ensures equilibrium is reached and samples configurations according to their equilibrium probabilities, adhering to detailed balance conditions. The likelihood of being in a particular state is determined by the system's energy and temperature.






