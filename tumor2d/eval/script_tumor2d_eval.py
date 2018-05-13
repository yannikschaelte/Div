import pyabc
import tumor2d
from pyabc.visualization import plot_kde_matrix
import matplotlib.pyplot as plt
import numpy as np
import math


limits = dict(log_division_rate=(-3, -1),
              log_division_depth=(1, 3),
              log_initial_spheroid_radius=(0, 1.2),
              log_initial_quiescent_cell_fraction=(-5, 0),
              log_ecm_production_rate=(-5, 0),
              log_ecm_degradation_rate=(-5, 0),
              log_ecm_division_threshold=(-5, 0))

prior = pyabc.Distribution(**{key: pyabc.RV('uniform', a, b - a)
                              for key, (a, b) in limits.items()})

h_loaded_adap = pyabc.History("sqlite:////media/sf_shared/tumor2d_2.db")
h_loaded_adap.id = 9

h_loaded_prev = pyabc.History("sqlite:////media/sf_shared/tumor2d.db")
h_loaded_prev.id = 1

h_loaded_prev1 = pyabc.History("sqlite:////media/sf_shared/tumor2d_1.db")
h_loaded_prev1.id = 1

h_loaded_dflt = pyabc.History("sqlite:///" + tumor2d.stored_data_db)

h_loaded_ctr = pyabc.History("sqlite:////media/sf_shared/tumor2d_adap.db")
# # kde plots

def kde_plots(name):
    
    if name == "dflt":
        h_loaded = h_loaded_dflt
    elif name == "prev":
        h_loaded = h_loaded_prev
    elif name == "prev1":
        h_loaded = h_loaded_prev1
    elif name == "adap":
        h_loaded = h_loaded_adap
    elif name == "ctr":
        h_loaded = h_loaded_ctr

    for t in range(h_loaded.max_t+1):
        df, w = h_loaded.get_distribution(m=0, t=t)
        plot_kde_matrix(df, w, limits=limits)
        plt.savefig(name + '_' + str(t))
        plt.close()
        print("done with " + name + "_" + str(t))

# # samples

pops_dflt = h_loaded_dflt.get_all_populations()
pops_prev = h_loaded_prev.get_all_populations()
pops_prev1 = h_loaded_prev1.get_all_populations()
pops_adap = h_loaded_adap.get_all_populations()
    
def samples():
    samples_dflt = pops_dflt['samples']
    samples_prev = pops_prev['samples']
    samples_prev1 = pops_prev1['samples']
    samples_adap = pops_adap['samples']
    
    t = np.arange(0,40)
    fig = plt.figure()
    plt.plot(t[:30],samples_dflt.values[:30], 'mo-', label="dflt")
    plt.plot(t[:len(samples_prev.values)],samples_prev.values, 'bo-', label="prev")
    plt.plot(t[:len(samples_prev1.values)],samples_prev1.values, 'go-', label="prev1")
    plt.plot(t[:len(samples_adap.values)],samples_adap.values, 'ro-', label="adap")
    plt.legend()
    fig.savefig("samples")
    plt.close()

# # epsilon

def epsilon():
    eps_dflt = pops_dflt['epsilon']
    eps_prev = pops_prev['epsilon']
    eps_prev1 = pops_prev1['epsilon']
    eps_adap = pops_adap['epsilon']

    t = np.arange(0,40)
    fig = plt.figure()
    plt.plot(t[:30],np.log(eps_dflt.values[:30]), 'mo-', label="dflt")
    plt.plot(t[:len(eps_prev.values)],np.log(eps_prev.values), 'bo-', label="prev")
    plt.plot(t[:len(eps_prev1.values)],np.log(eps_prev1.values), 'go-', label="prev1")
    plt.plot(t[:len(eps_adap.values)],np.log(eps_adap.values), 'ro-', label="adap")
    plt.legend()
    fig.savefig("eps")
    plt.close()

kde_plots("ctr")
#samples()
#epsilon()
