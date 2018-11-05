from model import *

simulate(model_need_bias, y_obs_error, "nada", "needBias_yObsError_nada")
simulate(model_need_bias, y_obs_error, "std", "needBias_yObsError_std")
simulate(model_need_bias, y_obs_error, "rmsd", "needBias_yObsError_rmsd")
