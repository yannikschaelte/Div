import numpy as np
import scipy as sp
import pyabc
import time


n_t = 100
n_r = 100

def model(p):
    time.sleep(1)
    return {'y' + str(j) : p['mean'] + np.random.randn(n_t) for j in range(0, n_r)}


limits = {'mean': (0, 1)}
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', bounds[0], bounds[1])
                              for key, bounds in limits.items()})

p_true = {'mean': 0.5}
y_obs = {'y' + str(j) : p_true['mean'] + np.zeros(n_t) for j in range(0, n_r)}


def distance(x, y):
    return np.power(x['y0'] - y['y0'], 2).sum()


pop_size = 1000
redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-lisa", port=8775)
multicore_sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=40)
