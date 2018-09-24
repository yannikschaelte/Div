import numpy as np
import scipy as sp
import scipy.integrate as integrate
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
    def uniform_dty(m):
        if m > prior_ub or m < prior_lb:
            return 0
        return 1 / (prior_ub - prior_lb)
    prior_val = uniform_dty(m)

    def normal_dty(m, m_true):
        return np.exp( - (m - m_true)**2 / (2 * std**2 / n_r ) )
    likelihood_val = normal_dty(m, th_true['m'])

    return likelihood_val * prior_val


# pyabc stuff
distance = pyabc.PNormDistance(p=2)
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=6)
max_nr_populations = 10
pop_size = 500

# visualize

def visualize(label, history, show_true=True, legend=True):
    t = history.max_t
    lb = -0.5
    ub = 0.5
    df, w = history.get_distribution(m=0, t=t)
    _, ax = plt.subplots()
    
    if show_true:
        integral = integrate.quad(pdf_true, prior_lb, prior_ub)[0]
        def pdf(x):
            return pdf_true(x) / integral

        xs = np.linspace(lb, ub, 10000)
        ys = []
        for x in xs:
            ys.append(pdf(x))
        ax.plot(xs, ys, '-', color='k', alpha=0.75)

    ax = pyabc.visualization.plot_kde_1d(df, w, xmin=lb, xmax=ub,
        x='m', numx=200, refval=th_true, ax=ax)
    
    if legend:
        if show_true:
            ax.legend(['true posterior', 'ABC posterior', 'true parameter'])
        else:
            ax.legend(['ABC posterior', 'true parameter'])

    plt.savefig(label + "_kde_1d_" + str(t) + ".png")
    plt.close()
