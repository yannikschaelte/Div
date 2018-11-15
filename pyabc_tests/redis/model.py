import numpy as np
import scipy as sp
import pyabc


n_t = 10000

def model(p):
    return {'y': p['mean'] + np.random.randn(n_t)}


limits = {'mean': (0, 1)}
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', bounds[0], bounds[1])
                              for key, bounds in limits.items()})

p_true = {'mean': 0.5}
y_obs = {'y': p_true['mean'] + np.zeros(n_t)}


def distance(x, y):
    return np.power(x['y'] - y['y'], 2).sum()


pop_size = 1000
redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host="icb-lisa", port=8775)
multicore_sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=20)
