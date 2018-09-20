import numpy as np
import scipy as sp
import pickle
import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import scipy.stats as stats
import os


# noise
std = 1
# number of replicates
n_r = 1000


def model(th):
    return th['m'] * np.ones(n_r) + std * np.random.randn(n_r)


def sumstat_all(y):
    return {'s' + str(j) : val for j, val in enumerate(y)}


def sumstat_samplemean(y):
    return {'mean': np.mean(y)}


# true parameter
th_true = {'m': 0}
# data
y_obs = np.zeros(n_r)
# sumstat
sumstat_all_obs = sumstat_all(y_obs)
sumstat_samplemean_obs = sumstat_samplemean(y_obs)
# prior
prior_lb = -4
prior_ub = 4
prior = pyabc.Distribution(
    **{'m': pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)})

# true pdf

def pdf_true(m):
    prior_val = prior.pdf(m)
    def normal_dty(m, m_true):
        return (1 / np.sqrt(2 * np.pi * std**2)
               * np.exp( - ((y_obs['m'] - y['m']) / noise) ** 2 / 2))
    likelihood_val = normal_dty(m, th_true['m'])

    return likelihood_val * prior_val


# pyabc stuff
distance = pyabc.PNormDistance(p=2)
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=6)
max_nr_populations = 1
pop_size = 500

# visualize

def visualize(label, history, show_true=True):
    t = history.max_t
    df, w = history.get_distribution(m=0, t=t)
    ax = pyabc.visualization.plot_kde_1d(df, w, xmin=prior_lb, xmax=prior_ub,
        x='m')
    plt.savefig(label + "_kde_1d_" + str(t) + ".png")
    plt.close()
