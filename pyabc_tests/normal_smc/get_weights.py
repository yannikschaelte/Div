import pyabc
import numpy as np

def get_weights():
    db_file = "sqlite:///db1.db"
    h = pyabc.History(db_file)
    h.id = 1
    n_effs = []
    for t in range(0, h.max_t + 1):
        d_weighted = h.get_weighted_distances(t)
        w = np.array(d_weighted['w'])
        n_eff = np.sum(w)**2 / np.sum(w**2)
        n_effs.append(n_eff)
        print(n_eff)
    return n_effs
