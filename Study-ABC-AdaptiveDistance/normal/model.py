import scipy as sp
import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import os
import tempfile
import logging
import sys
import copy
import numpy as np


# for debugging
logger = logging.getLogger('DistanceFunction')
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)

fh = logging.FileHandler("logging.log")
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)


# model definition
def model_need_ad(p):
    return {'s0': p['p0'] + 0.01 * sp.randn(),
            's1': p['p1'] + 1 * sp.randn(),
            's2': 2 + 100 * sp.randn()}

def model_need_bias(p):
    s = model_need_ad(p)
    s['s2'] = 2 + 0.0001 * sp.randn()

# true model parameter
p_true = {'p0': 4,
          'p1': 5}

# observed summary statistics
y_obs = {'s0': p_true['p0'],
         's1': p_true['p1'], 
         's2': 2}

y_obs_error = copy.deepcopy(y_obs)
y_obs_error['s2'] = 8

# prior distribution
limits = {'p0': (0, 10), 'p1': (0, 10)}
prior = pyabc.Distribution(**{key:pyabc.RV('uniform', bounds[0], bounds[1])
                              for key, bounds in limits.items()})

# further variables
max_nr_populations = 8
sampler = pyabc.sampler.MulticoreEvalParallelSampler(n_procs=16)

# visualization
def visualize(label, history):
    fig, ax = plt.subplots()
    for t in range(1, history.max_t + 1):
        df, w = history.get_distribution(m=0, t=t)
        pyabc.visualization.plot_kde_matrix(df, w, limits=limits, 
                                            refval = p_true)
        plt.savefig("kde_" + label + "_" + str(t) + ".png")
        plt.close()
    
    samples = history.get_all_populations()['samples']
    fig = plt.figure()
    plt.plot(samples.values, 'o-', label="samples")
    plt.savefig("samples_" + label + "_" + str(t) + ".png")
    plt.close()


def viz_eps(list_h, list_label):
    list_eps = []
    for h in list_h:
        list_eps.append(np.array(h.get_all_populations()['epsilon']))
    _, ax = plt.subplots()
    for ix, eps_schedule in enumerate(list_eps):
        ax.plot(np.log(eps_schedule[1:]), 'x-', label=list_label[ix])
    plt.xlabel("Population index")
    plt.ylabel("Log(Epsilon)")
    plt.legend()
    plt.savefig("viz_eps.png")


def viz_samples(list_h, list_label):
    list_samples = []
    for h in list_h:
        list_samples.append(np.array(h.get_all_populations()['samples']))
    _, ax = plt.subplots()
    for ix, sample_schedule in enumerate(list_samples):
        ax.plot(np.log(sample_schedule[1:]), 'x-', label=list_label[ix])
    plt.xlabel("Population index")
    plt.ylabel("Log(#Samples)")
    plt.legend()
    plt.savefig("viz_samples.png")


# simulation function

def simulate(model, y_obs, method, label):

    if method == "nada":
        distance = pyabc.PNormDistance(p=2)
    else:
        if method == "mad":
            scale_function = pyabc.distance_functions.median_absolute_deviation
        elif method == "mado":
            scale_function = pyabc.distance_functions.median_absolute_deviation_to_observation
        elif method == "std":
            scale_function = pyabc.distance_functions.standard_deviation
        elif method == "bias":
            scale_function = pyabc.distance_functions.bias
        elif method == "rmsd":
            scale_function = pyabc.distance_functions.root_mean_square_deviation
        elif method == "cmad":
            scale_function = pyabc.distance_functions.combined_median_absolute_deviation
        distance = pyabc.AdaptivePNormDistance(
                p=2, adaptive=True, scale_function=scale_function)

    abc = pyabc.ABCSMC(model, prior, distance, sampler=sampler)
    db_path = "sqlite:///db_" + label + ".db"
    abc.new(db_path, y_obs)
    h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
    visualize(label, h)
