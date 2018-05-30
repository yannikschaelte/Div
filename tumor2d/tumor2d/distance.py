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


class ReweightedTumor2DDistance(Tumor2DDistance):
    
    def __init__(self, variances: dict):
        super().__init__(variances)
        self.inv_variances['growth_curve'] /= 20
        self.inv_variances['proliferation_profile'] /= 1000
        self.inv_variances['extra_cellular_matrix_profile'] /= 1000


class AdaptiveTumor2DDistance(AdaptivePNormDistance):
    

    def __init__(self, adaptive=True):
        super().__init__(p=2, 
                         adaptive=adaptive,
                         scale_type=AdaptivePNormDistance.SCALE_TYPE_C_MAD)
    
    def initialize(self, t, sample_from_prior, x_0):
        sum_stats = []
        for sum_stat in sample_from_prior:
            sum_stats.append(normalize_sum_stats(sum_stat))
        x_0 = normalize_sum_stats(x_0)
        super().initialize(t, sum_stats, x_0)

    def update(self, t, all_sum_stats, x_0):
        sum_stats = []
        for sum_stat in all_sum_stats:
            sum_stats.append(normalize_sum_stats(sum_stat))
        x_0 = normalize_sum_stats(x_0)
        super().update(t, sum_stats, x_0)

    def __call__(self, t, x, y):
        x = normalize_sum_stats(x)
        y = normalize_sum_stats(y)
        return super().__call__(t, x, y)


class ReweightedAdaptiveTumor2DDistance(AdaptiveTumor2DDistance):

    def update(self, t, all_sum_stats, x_0):
        super().update(t, all_sum_stats,, x_0)
        for key in self.w[t]:
            a, b = key
            if a == 'growth_curve':
                self.w[t][key] /= 20
            elif a == 'proliferation_profile' or a == 'extra_cellular_matrix_profile':
                self.w[t][key] /= 1000


def normalize_sum_stats(x):
    x_flat = {}
    for key, value in x.items():
        for j in range(len(value)):
            x_flat[(key, j)] = value[j]
    return x_flat
