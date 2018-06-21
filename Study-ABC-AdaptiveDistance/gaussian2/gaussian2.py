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
    return {'s0': p['th0'] + 1 + 0.1*sp.randn(),
            's1': p['th1'] + 2 + 2*sp.randn()}


# true model parameter
th_true = {'th0': 2, 'th1': 5}

# observed summary statistics
observation = {'s0': th_true['th0'] + 1, 's1': th_true['th1'] + 2}

# prior distribution
prior = pyabc.Distribution(th0=pyabc.RV('uniform', 0, 10),
                           th1=pyabc.RV('uniform', 0, 10))

# further variables
max_nr_populations = 10


# visualization
def visualize(label, history):
    fig, ax = plt.subplots()
    t = history.max_t
    df, w = history.get_distribution(m=0, t=t)
    pyabc.visualization.plot_kde_matrix(df, w, limits={key: (0, 10) for key in ['th0', 'th1']}, 
                                    refval = th_true)
    plt.savefig("plt_kde_matrix_" + label + "_" + str(t))

    samples = history.get_all_populations()['samples']
    fig = plt.figure()
    plt.plot(samples.values, 'o-', label="samples")
    plt.savefig("samples_" + label + "_" + str(t))


def simulate(method):
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

    abc = pyabc.ABCSMC(model, prior, distance)
    db_path = "sqlite:///db.db"
    abc.new(db_path, observation)
    h = abc.run(minimum_epsilon=0, max_nr_populations=max_nr_populations)
    visualize(method, h)

# command line arguments:
# nada, mad, mado, std, bias, rmsd, cmad
method = sys.argv[1]
simulate(method)
