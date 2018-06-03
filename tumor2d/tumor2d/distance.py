from .simulate import nr_valid
import numpy as np
import pyabc
import logging
df_logger = logging.getLogger("DistanceFunction")

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


class AdaptiveTumor2DDistance(pyabc.AdaptivePNormDistance):
    

    def __init__(self, adaptive=True):
        super().__init__(p=2, 
                         adaptive=adaptive,
                         scale_function=pyabc.distance_functions.median_absolute_deviation)
    
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
        super().update(t, all_sum_stats, x_0)
        for key in self.w[t]:
            a, b = key
            if a == 'growth_curve':
                self.w[t][key] /= 20
            elif a == 'proliferation_profile' or a == 'extra_cellular_matrix_profile':
                self.w[t][key] /= 1000


class LessWeightsAdaptiveTumor2DDistance(pyabc.DistanceFunction):
    """
    Only one weight for each of the data types growth_curve,
    proliferation_profile, extra_cellular_matrix_profile.
    """
    
    def __init__(self):
        self.require_initialize = True
        self.w = {}

    def initialize(self, t, sample_from_prior, x_0):
        self._update(t, sample_from_prior, x_0)

    def update(self, t, all_sum_stats, x_0):
        self._update(t, sample_from_prior, x_0)

    def _update(self, t, sample_from_prior, x_0):
        scales = {}
        scale_fun = pyabc.median_absolute_deviation
        n = len(sample_from_prior)
        w = {}
        for key in sample_from_prior[0]:
            max_len = max(len(sample_from_prior[j][key]) for j in range(n))
            for j in range(max_len):
                current_list = []
                for ss in sample_from_prior:
                    if len(ss[key]) > j:
                        current_list.append(ss[key][j])
                scale = scale_fun(current_list)
                scales.setdefault(key, []).append(scale)
            scale = np.mean(scales)
            if np.isclose(scale, 0):
                w[key] = 0
            else:
                w[key] = 1 / scale
        mean_weight = np.mean(list(w.values()))
        for key in w:
            w[key] /= mean_weight

        self.w[t] = w
        df_logger.debug("update distance weights = {}".format(self.w[t]))

    def configure_sampler(self, sampler):
        sampler.sample_factory.record_all_sum_stats = True

    def get_config(self):
        return {"name": self.__class__.__name__}

    def __call__(self, t, x, y):
        w = self.w[t]
        d = sum(
                sum(
                    pow(abs(w[key]*(x[key][j]-y[key][j])), 2) 
                for j in range(min(len(x[key]), len(y[key]))))
            for key in w)

        return d

def normalize_sum_stats(x):
    x_flat = {}
    for key, value in x.items():
        for j in range(len(value)):
            x_flat[(key, j)] = value[j]
    return x_flat
