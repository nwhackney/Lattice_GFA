Markov chain Monte Carlo code to study finite-temperature behavior of lattice model of geometrically frustrated assembly.
Details of model and Monte Carlo algorithm outlined by Hackney et. al. PRX 2023 (https://journals.aps.org/prx/pdf/10.1103/PhysRevX.13.041010)

Example input script "params.json" included and can be interpreted as follows:

"mc_params" specify details of Monte Carlo procedure.
  Number of MC sweeps to make before outputting data set by "step_size" and total number of these blocks to run set by "num_steps."
  "seed" is for random number generation
  "restart:" specifies if initial config should be read from checkpoint file (i.e. "restart: 1") or if new random state should be generated (i.e. "restart:0")

"xy_params" specify parameters for the model outlined in paper referenced above.

  beta defines the thermal parameter beta J
  sigma defines the cohesion-to-stiffness parameter Sigma/J
  f defines the frustration and ranges from 0 (unfrustrated) to 1/2 (fully frustrated)
  length defines the length, L, of the L by L square lattice
  and Phi defines the fraction of occupied lattice sites (i.e. number of particles / L*L)
