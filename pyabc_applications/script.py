import numpy as np
from concurrent.futures import ThreadPoolExecutor
from pyabc import Distribution, RV, ABCSMC
from pyabc.sampler import ConcurrentFutureSampler
from tumor2d import simulate, log_model, load_default, Tumor2DDistance


# MODEL

model = log_model

# PRIOR

limits = dict(log_division_rate=(-3, -1),
              log_division_depth=(1, 3),
              log_initial_spheroid_radius=(0, 1.2),
              log_initial_quiescent_cell_fraction=(-5, 0),
              log_ecm_production_rate=(-5, 0),
              log_ecm_degradation_rate=(-5, 0),
              log_ecm_division_threshold=(-5, 0))

prior = Distribution(**{key: RV("uniform", a, b - a)
                        for key, (a,b) in limits.items()})

# OBSERVATIONS

# observation = simulate(division_rate=4.17e-2,
#                        initial_spheroid_radius=1.2e1,
#                        initial_quiescent_cell_fraction=7.5e-1,
#                        division_depth=100,
#                        ecm_production_rate=5e-3,
#                        ecm_degradation_rate=8e-4,
#                        ecm_division_threshold=1e-2)

# as default observation we consider the mean of 100 samples generated with the
# above reference parameters

data = load_default()
data_mean = data[1]
data_var = data[2]

observation = data_mean

# DISTANCE

distance = Tumor2DDistance(data_var)

# SAMPLER

pool = ThreadPoolExecutor(max_workers=2)
sampler = ConcurrentFutureSampler(pool)

# ABC

abc = ABCSMC(models=model,
             parameter_priors=prior,
             distance_function=distance,
             population_size=3,
             sampler=sampler)

abc.new(db="sqlite:////tmp/test.db",
        observed_sum_stat=observation)

history = abc.run(max_nr_populations=1,
                  minimum_epsilon=0)
