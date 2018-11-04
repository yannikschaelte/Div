import pyabc
from model import viz_eps, viz_samples

list_label = ["Nothing", "STD", "RMSD"]
name_base = "sqlite:///db_needAd_yObs_"
list_h = [pyabc.History(name_base + "nada.db"),
          pyabc.History(name_base + "std.db"),
          pyabc.History(name_base + "rmsd.db")]

viz_eps(list_h, list_label)
viz_samples(list_h, list_label)
