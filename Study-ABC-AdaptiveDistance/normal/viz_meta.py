import pyabc
from model import viz_eps, viz_samples

list_label = ["Nothing", "STD", "RMSD"]
name_base = "sqlite:///db_needAd_yObs_"
list_h = [pyabc.History(name_base + "nada.db"),
          pyabc.History(name_base + "std.db"),
          pyabc.History(name_base + "rmsd.db")]
label="viz_needAd_yObs_"
viz_eps(label, list_h, list_label)
viz_samples(label, list_h, list_label)


list_label = ["Nothing", "STD", "RMSD"]
name_base = "sqlite:///db_needAd_yObsError_"
list_h = [pyabc.History(name_base + "nada.db"),
          pyabc.History(name_base + "std.db"),
          pyabc.History(name_base + "rmsd.db")]
label="viz_needAd_yObsError_"
viz_eps(label, list_h, list_label)
viz_samples(label, list_h, list_label)


list_label = ["Nothing", "STD", "RMSD"]
name_base = "sqlite:///db_needBias_yObsError_"
list_h = [pyabc.History(name_base + "nada.db"),
          pyabc.History(name_base + "std.db"),
          pyabc.History(name_base + "rmsd.db")]
label="viz_needBias_yObsError_"
viz_eps(label, list_h, list_label)
viz_samples(label, list_h, list_label)
