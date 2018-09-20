import numpy as np
import pyabc.visualization
import os
import tempfile
import matplotlib.pyplot as plt
import pickle

# model


def gk(p):
    # extract parameters
    a = p['A']
    b = p['B']
    g = p['g']
    k = p['k']
    c = 0.8

    # sample from normal distribution
    z = np.random.normal(0, 1)

    # to sample from gk distribution
    return _z_to_gk(a, b, g, k, c, z)


def _z_to_gk(a, b, g, k, c, z):
    e = np.exp(-g*z)
    return a + b * (1 + c * (1-e) / (1+e)) * (1 + z**2)**k * z


# number of replicates
n_r = 1000
octile = n_r // 8
order_statistics_indices = [octile * j for j in range(1, 8)]


def model(p):
    y = [gk(p) for _ in range(0, n_r)]
    return y


# sum stats


def _ordstat(y: list):
    y = sorted(y)
    ordstat = {j: y[index] 
               for j, index in enumerate(order_statistics_indices)}
    return ordstat


def sumstat_all(y):
    sumstat = {'y' + str(j): val for j, val in enumerate(y)}
    return sumstat


def sumstat_dp(y):
    """
    According to Drovandi-Pettitt.
    """
    ordstat = _ordstat(y)
    sumstat = {}
    sumstat['A'] = ordstat[3]
    sumstat['B'] = ordstat[5] - ordstat[1]
    sumstat['g'] = (ordstat[5] + ordstat[1] - 2 * ordstat[3]) / sumstat['B']
    sumstat['k'] = (ordstat[6] - ordstat[4] + ordstat[2] - ordstat[0]) \
        / sumstat['B']
    return sumstat


# true parameters
p_true = {'A': 3.0,
          'B': 1.0,
          'g': 2.0,
          'k': 0.5}

# observed data
_y_obs = None

def get_y_obs():
    global _y_obs
    if _y_obs is not None:
        return _y_obs

    y_file = "y.dat"
    try:
        y_obs = pickle.load(open(y_file, 'rb'))
    except Exception:
        y_obs = model(p_true)
        #y_obs_arr = []
        #for _ in range(0, 1000):
        #    y_obs_arr.append(model(p_true))
        #y_obs_arr = np.array(y_obs_arr)
        #y_obs = np.mean(y_obs_arr, axis=0)
        pickle.dump(y_obs, open(y_file, 'wb'))
    
    _y_obs = y_obs
    return y_obs


y_obs = get_y_obs()


# observed summary statistics
sumstat_all_obs = sumstat_all(y_obs)
sumstat_dp_obs = sumstat_dp(y_obs)

# prior
prior_lb = 0
prior_ub = 5
prior = pyabc.Distribution(**{key: pyabc.RV('uniform', prior_lb, prior_ub)
                              for key in p_true})

# pyabc stuff
distance = pyabc.AdaptivePNormDistance()
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=6)
max_nr_populations = 10
pop_size = 100

# visualize

def visualize(label, history):
    t = history.max_t
    df, w = history.get_distribution(m=0, t=t)
    ax = pyabc.visualization.plot_kde_matrix(
            df, w,
            limits={key: (prior_lb, prior_ub)
                    for key in p_true},
            refval=p_true)
    plt.suptitle("Posterior KDE plot")
    plt.savefig(label + "_kde_matrix_" + str(t) + ".png")
    plt.close()
