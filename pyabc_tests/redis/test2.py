import numpy as np
from model import *
import pyabc
from pyabc.external import R


r = R("r_model.r")
model = r.model("model")
sumstat = r.summary_statistics("sumstat")
distance = r.distance("distance")
y_obs = r.observation("y_obs")


abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=sumstat,
                   sampler=redis_sampler,
                   population_size=pop_size)


db = "sqlite:///db1.db"
abc.new(db, observed_sum_stat=y_obs)
abc.run(minimum_epsilon=0.0, max_nr_populations=30, min_acceptance_rate=0.001)
