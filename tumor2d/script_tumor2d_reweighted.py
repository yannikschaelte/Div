import numpy as np
import time
from pyabc import Distribution, RV, ABCSMC, LocalTransition, MultivariateNormalTransition
from pyabc.sampler import *
from pyabc.populationstrategy import *
from tumor2d import simulate, log_model, load_default, Tumor2DDistance, ReweightedTumor2DDistance

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

distance = ReweightedTumor2DDistance(data_var)
#distance = AdaptivePNormDistance(p=2, 
#								 use_all_w=True, 
#								 adaptive=True, 
#								 scale_type=AdaptivePNormDistance.SCALE_TYPE_SD)

# SAMPLER

sampler = RedisEvalParallelSampler(host="wastl", port=8765)
# sampler = SingleCoreSampler()

# POPULATION STRATEGY

# population_size = AdaptivePopulationSize(start_nr_particles=500, mean_cv=0.25, n_bootstrap=5)
population_size = 300

# TRANSITION STRATEGY

# transition = LocalTransition(k=100)
transition = MultivariateNormalTransition()

# PREPARE ABC

# for real run: 100 cores, initial population size 500

abc = ABCSMC(models=model,
             parameter_priors=prior,
             distance_function=distance,
             population_size=population_size,
             transitions=transition,
             sampler=sampler)

db_file = 'sqlite:////home/icb/yannik.schaelte/Div/tumor2d/tumor2d_reweighted.db'

abc.new(db=db_file, observed_sum_stat=observation)

# RUN ABC

start_time = time.time()
print("start abc.run, ", time.asctime(time.localtime(start_time)))

history = abc.run(max_nr_populations=30, minimum_epsilon=0)

print(f"done abc.run, {time.time() - start_time:.2f}s")
