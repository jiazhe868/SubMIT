Subevent Inversion Code

Overview
This code implements a multiple-subevent inversion method to constrain the source parameters of a large earthquake. The number of subevents is determined iteratively based on improvements in data fit. Each subevent is characterized by five nonlinear parameters (horizontal location, depth, centroid time, and centroid duration) and five linear parameters (deviatoric moment tensor components). The first subevent’s location is fixed at the hypocenter to prevent arbitrary shifts.

The inversion employs a two-stage algorithm:
Nonlinear parameter search using Metropolis-Hastings Markov Chain Monte Carlo (MCMC).
Linear inversion for subevent moment tensors for each set of nonlinear parameters.
A Bayesian framework is used to incorporate data errors and model priors, ensuring uncertainty quantification based on posterior probability distributions.

Requirements
High-performance computing resources (multi-core processing).
Par.file and search_par.file contain key inversion parameters.

Key Control Files
Par.file (Controls the inversion process)
Defines time windows, sampling rates, station data, and inversion weights.
Specifies frequency ranges, regularization parameters, and moment tensor scaling.
Includes source model constraints (depth, rupture velocity, etc.).

Input.model (Initial Subevent Model)
Each row represents a subevent with the following columns:
Centroid time (s)
EW location (km)
NS location (km)
Duration (s)
Rupture velocity (km/s) (0.01 ≈ point source, but avoid using 0)
Rupture directivity (°) (not meaningful for point sources)
Centroid depth (km)

search_par.file (Controls MCMC Search)
Burn_in and nsample: Typically set to n_subevent * 1000 but may need longer runs (10,000 burn-in, 5,000 samples for complex events).
Search boundaries: Defined for location, depth, duration, and centroid time, based on aftershock distribution.
Number of subevents (neq_min and neq_max): Should be set identically here.
