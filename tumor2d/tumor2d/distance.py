from .simulate import nr_valid
import numpy as np
from pyabc import AdaptivePNormDistance

class Tumor2DDistance:
    __name__ 
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


def non_weighted_tumor2ddistance(x: dict, y: dict) -> float:
    length = {key: min([len(x[key]),
                        len(y[key])]) for key in y}
    return sum(np.sum((x[key][:length[key]] - y[key][:length[key]])**2) for key in length.keys())


class AdaptiveTumor2DDistance(AdaptivePNormDistance):
    

    def __init__(self):
        super().__init__(p=2, 
                         use_all_w=True, 
                         adaptive=True, 
                         scale_type=AdaptivePNormDistance.SCALE_TYPE_SD)
    
    def initialize(self, t, sample_from_prior):
        sum_stats = []
        for sum_stat in sample_from_prior:
            sum_stats.append(normalize_sum_stats(sum_stat))
        super().initialize(t, sum_stats)

    def update(self, t, all_sum_stats):
        sum_stats = []
        for sum_stat in all_sum_stats:
            sum_stats.append(normalize_sum_stats(sum_stat))
        super().update(t, sum_stats)

    def __call__(self, t, x, y):
        x = normalize_sum_stats(x)
        y = normalize_sum_stats(y)
        return super().__call__(t, x, y)


def normalize_sum_stats(x):
    x_flat = {}
    for key, value in x.items():
        for j in range(len(value)):
            x_flat[(key, j)] = value[j]
    return x_flat
