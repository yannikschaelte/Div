from model import *
import pyabc

db_file = "sqlite:///db1.db"
abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=sumstat_all,
                   population_size=pop_size,
                   sampler=sampler)
#abc.new(db=db_file, observed_sum_stat=sumstat_all_obs)
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
h = pyabc.History(db_file)
visualize("pic1", h)
