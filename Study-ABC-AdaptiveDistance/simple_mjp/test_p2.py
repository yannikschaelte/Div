from model import *
import pyabc

db_file= "sqlite:///db2.db"

abc = pyabc.ABCSMC(models=model_random,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=sumstat_p,
                   population_size=pop_size,
                   sampler=sampler)
#abc.new(db=db_file, observed_sum_stat=sumstat_p(y_obs))
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
h = pyabc.History(db_file)
#h.id = 1
visualize("viz2_" + str(n_r), h)
