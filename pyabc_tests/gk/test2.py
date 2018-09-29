from model import *
import pyabc
import pyabc.visualization


db_file = "sqlite:///db2.db"
abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=sumstat_dp,
                   population_size=pop_size,
                   sampler=sampler)
#abc.new(db=db_file, observed_sum_stat=sumstat_dp_obs)
#h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
h = pyabc.History(db_file)

visualize("pic2", h)
