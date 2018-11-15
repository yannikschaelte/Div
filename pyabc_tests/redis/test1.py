import numpy as np
from model import *
import pyabc


abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=lambda x: x,
                   sampler=multicore_sampler,
                   population_size=pop_size)


db = "sqlite:///db1.db"
abc.new(db, observed_sum_stat=y_obs)
abc.run(minimum_epsilon=0.0, max_nr_populations=30, min_acceptance_rate=0.01)
