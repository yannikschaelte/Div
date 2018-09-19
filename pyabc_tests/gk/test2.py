from model import *
import pyabc
import pyabc.visualization


db_file = "sqlite:///db2.db"
abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=distance,
                   summary_statistics=sumstat_dp,
                   population_size=500)
abc.new(db=db_file, observed_sum_stat=sumstat_all_obs)
h = abc.run(minimum_epsilon=0, max_nr_populations=10)

visualize("test2", h)
