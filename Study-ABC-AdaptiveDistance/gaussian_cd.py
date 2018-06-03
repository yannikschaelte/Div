import scipy as sp
import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt
import os
import tempfile
import logging
import sys


# for debugging
df_logger = logging.getLogger('DistanceFunction')
df_logger.setLevel(logging.DEBUG)


# model definition
def model(p):
    return {'ss1': p['theta'] + 1 + 0.1*sp.randn(),
            'ss2': 2 + 2*sp.randn()}


# true model parameter
theta_true = 3

# observed summary statistics
observation = {'ss1': theta_true + 1, 'ss2': 8}

# prior distribution
prior = pyabc.Distribution(theta=pyabc.RV('uniform', 0, 10))

# further variables
max_nr_populations = 8

# visualization
def visualize(label, history):
    fig, ax = plt.subplots()
    t = history.max_t
    df, w = history.get_distribution(m=0, t=t)
    pyabc.visualization.plot_kde_1d(df, w, xmin=0, xmax=10,
                                    x='theta', ax=ax, 
                                    refval = {'theta': theta_true},
                                    label="PDF t={}".format(t))
    plt.savefig("plt_kde_1d_" +label + "_" + str(t))

    samples = history.get_all_populations()['samples']
    fig = plt.figure()
    plt.plot(samples.values, 'o-', label="samples")
    plt.savefig("samples_" + label + "_" + str(t))


def simulate(method):
    if method == "nada":
        distance = pyabc.PNormDistance(p=2)
    elif method == "mad":
        distance = pyabc.AdaptivePNormDistance(
                p=2,
                adaptive=True,
                scale_function = pyabc.distance_functions.median_absolute_deviation)
    elif method == "cmad":
        distance = pyabc.AdaptivePNormDistance(
                p=2,
                adaptive=True,
                scale_function = pyabc.distance_functions.centered_median_absolute_deviation)

    abc = pyabc.ABCSMC(model, prior, distance)
    db_path = "sqlite:///db.db"
    abc.new(db_path, observation)
    h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
    visualize(method, h)

# command line arguments:
# nada, mad, cmad
method = sys.argv[1]
simulate(method)
