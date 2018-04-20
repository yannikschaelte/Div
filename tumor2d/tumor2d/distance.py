from .simulate import nr_valid
import numpy as np
from pyabc import AdaptivePNormDistance


class Tumor2DDistance:
    """
    No subclass of pyabc.Distance, so will be converted to SimpleFunctionDistance.
    """
    def __init__(self, variances: dict):
        self.variances = {key: val[:nr_valid(val)]
                          for key, val in variances.items()}
        self.inv_variances = {}
        for key, val in self.variances.items():
            inv = np.zeros(len(val))
            inv[val != 0] = 1 / val[val != 0]
            self.inv_variances[key] = inv
    
    def initialize(self, sample_from_prior):
        pass
    
    def __call__(self, x: dict, y: dict) -> float:
        length = {key: min([len(x[key]),
                            len(self.variances[key]),
                            len(y[key])]) for key in y}
        return sum(np.sum((x[key][:length[key]] - y[key][:length[key]])**2
                              * self.inv_variances[key][:length[key]])
                   for key in length.keys())
    
    def get_config(self):
        return {}


class AdaptiveTumor2DDistance(AdaptivePNormDistance):
    """
    Make it adaptive.
    """

    def initialize(self, t, sample_from_prior):
        normalized_sum_stats = []
        for sum_stats in sample_from_prior:
            normalized_sum_stats.append(normalize_sum_stats(sum_stats))
        super().initialize(t, normalized_sum_stats)

    def update(self, t, all_sum_stats):
        normalized_sum_stats = []
        for sum_stats in all_sum_stats:
            normalized_sum_stats.append(normalize_sum_stats(sum_stats))
        super().update(t, all_sum_stats)

    def __call__(self, t, x, y):
        x = normalize_sum_stats(x)
        y = normalize_sum_stats(y)
        super().call(t, x, y)


def normalize_sum_stats(x):
    """
    Flatten.
    """
    return x