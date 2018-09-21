import numpy as np
import scipy as sp
import scipy.integrate as integrate
import pickle
import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import scipy.stats as stats
import os


def mean(p):
    return 2 * (p - 2) * p * (p + 2)


def std(p):
    return np.sqrt(0.1 + p**2)


def model(th):
    p = th['th0']
    y = mean(p) + std(p) * np.random.randn()
    return y


def sumstat(y):
    s = {'s' : y}
    return s


# data
y_obs = 2
# sumstat
sumstat_obs = sumstat(y_obs)
# prior
prior_lb = -5
prior_ub = 5
prior = pyabc.Distribution(
    **{'th0': pyabc.RV('uniform', prior_lb, prior_ub - prior_lb)})

# true pdf

def pdf_true(p):
    def uniform_dty(p):
        if p > prior_ub or p < prior_lb:
            return 0
        return 1 / (prior_ub - prior_lb)
    prior_val = uniform_dty(p)

    def normal_dty(y_obs, mean, std):
        return np.exp( - (y_obs - mean)**2 / (2 * std**2) 
    likelihood_val = normal_dty(y_obs, mean(p), std(p))

    return likelihood_val * prior_val


# pyabc stuff
distance = pyabc.PNormDistance(p=1)
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=20)
max_nr_populations = 20
pop_size = 500

# visualize

def visualize(label, history, show_true=True):
    fig, ax = plt.subplots()
    for t in range(0, history.max_t + 1):
        df, w = history.get_distribution(m=0, t=t)
        ax = pyabc.visualization.plot_kde_1d(
            df, w, xmin=prior_lb, xmax=prior_ub,
            x='th0', numx=200, label="PDF t={}".format(t), 
            ax=ax)

    if show_true:
        integral = integrate.quad(pdf_true, prior_lb, prior_ub)[0]
        
        def pdf(x):
            return pdf_true(x) / integral

        xs = np.linspace(prior_lb, prior_ub, 200)
        ys = []
        for x in xs:
            ys.append(pdf(x))

        integral2 = integrate.quad(pdf, prior_lb, prior_ub)[0]
        print(integral, integral2)
        ax.plot(xs, ys, '-', color='0.75', label="True PDF")
    ax.legend()
    plt.savefig(label + "_kde_1d_" + str(t) + ".png")
    plt.close()
