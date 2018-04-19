from .simulate import nr_valid
import numpy as np
from pyabc import AdaptivePNormDistance

class Tumor2DDistance:
	"""
	No subclass of pyabc.Distance, so will be converted to SimpleFunctionDistance.
	"""
    __name__ 
    def __init__(self, variances: dict):
        self.variances = {key: val[:nr_valid(val)]
                          for key, val in variances.items()}
        self.inv_variances = {}
        for key, val in self.variances.items():
            inv = np.zeros(len(val))
            inv[val != 0] = 1 / val[val != 0]
            self.inv_variances[key] = inv
    
    def initialize(self, t, sample_from_prior):
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
		
		
class Tumor2DPNormDistance(DistanceFunction):
    """
    Use weighted p-norm

    .. math::

        d(x, y) =\
         \\left[\\sum_{i} \\left w_i| x_i-y_i \\right|^{p} \\right]^{1/p}

    to compute distances between sets of summary statistics. E.g. set p=2 to
    get a Euclidean distance.

    Parameters
    ----------

    p: float
        p for p-norm. Required p >= 1, p = math.inf allowed (infinity-norm).

    w: dict
        Weights. Dictionary indexed by time points. Each entry contains a
        dictionary of numeric weights, indexed by summary statistics labels.
        If none is passed, a weight of 1 is considered for every summary
        statistic. If no entry is available in w for a given time point,
        the maximum available time point is selected.

    use_all_w: bool
        True: When checking for acceptance, check against all distances
        encoded.
        False: Only check against the distance at the passed time point.
    """

    def __init__(self,
                 p: float=2,
                 w: dict=None,
                 use_all_w: bool=True):
        super().__init__(require_initialize=False)

        if p < 1:
            raise ValueError("It must be p >= 1")
        self.p = p
        self.w = w
        self.use_all_w = use_all_w

    def __call__(self,
                 t: int,
                 x: dict,
                 y: dict) -> float:
        # make sure weights are initialized
        if self.w is None:
            self._set_default_weights(t, x.keys())

        # select time point
        if t not in self.w.keys():
            t = max(self.w.keys())

        # extract for time point
        w = self.w[t]

        if self.p == math.inf:
            d = max(abs(w[key]*(x[key]-y[key]))
                    for key in w.keys())
        else:
            d = pow(
                sum(pow(abs(w[key]*(x[key]-y[key])), self.p)
                    for key in w.keys()),
                1/self.p)

        return d

    def accept(self,
               t: int,
               eps: Epsilon,
               x: dict,
               x_0: dict) -> (float, bool):
        d = self(t, x, x_0)
        accept = d <= eps(t)

        if accept and self.use_all_w:
            # also check against all previous distances and acceptance criteria
            for t_prev in self.w.keys():
                if t_prev < t:
                    d_prev = self(t_prev, x, x_0)
                    accept = d_prev <= eps(t_prev)
                    if not accept:
                        break

        return d, accept

    def _set_default_weights(self,
                             t: int,
                             summary_statistics_keys):
        """
        Init weights to 1 for every summary statistic.
        """
        self.w = {t: {k: 1 for k in summary_statistics_keys}}


class AdaptiveTumor2DPNormDistance(Tumor2DPNormDistance):
    """
    Use a weighted p-norm to compute distances between sets of summary
    statistics.

    Parameters
    ----------

    p: float
        p for p-norm. Required p >= 1, p = math.inf allowed (infinity-norm).

    adaptive: bool
        True: Adapt distance after each iteration.
        False: Adapt distance only once at the beginning in initialize().
            This corresponds to a pre-calibration.

    scale_type: int
        What measure to use for deviation. Values as in SCALE_... constants.
    """

    # mean absolute deviation
    SCALE_TYPE_MAD = 0

    # standard deviation
    SCALE_TYPE_SD = 1

    def __init__(self,
                 p: float=2,
                 use_all_w: bool=True,
                 adaptive: bool=True,
                 scale_type: int=SCALE_TYPE_MAD):
        # call p-norm constructor
        super().__init__(p=p,
                         w=None,
                         use_all_w=use_all_w)
        self.require_initialize = True
        self.adaptive = adaptive
        self.scale_type = scale_type

    def configure_sampler(self,
                          sampler: Sampler):
        """
        Make the sampler return also rejected summary statistics if required,
        because these are needed to get a better estimate of the summary
        statistic variabilities.

        Parameters
        ----------

        sampler: Sampler
            The sampler employed.
        """

        super().configure_sampler(sampler)
        if self.adaptive:
            sampler.sample_factory.record_all_sum_stats = True

    def initialize(self,
                   t: int,
                   sample_from_prior: List[dict]):
        """
        Initialize weights.
        """

        super().initialize(t, sample_from_prior)
        # update weights from samples
        self._update(t, sample_from_prior)

    def update(self,
               t: int,
               all_sum_stats: List[dict]):
        """
        Update weights based on all simulations.
        """

        if not self.adaptive:
            return False

        self._update(t, all_sum_stats)

        return True

    def _update(self,
                t: int,
                all_sum_stats: List[dict]):
        """
        Here the real update of weights happens.
        """

        # retrieve keys
        keys = all_sum_stats[0].keys()

        # make sure w_list is initialized
        if self.w is None:
            self.w = {}

        n = len(all_sum_stats)

        # to-be-filled-and-appended weights dictionary
        w = {}

        for key in keys:
			
            # prepare list for key
            current_list = []
            for j in range(n):
                current_list.append(all_sum_stats[j][key])

            # compute weighting
            if self.scale_type == AdaptivePNormDistance.SCALE_TYPE_MAD:
                val = median_absolute_deviation(current_list)
            elif self.scale_type == AdaptivePNormDistance.SCALE_TYPE_SD:
                val = standard_deviation(current_list)
            else:
                raise Exception(
                    "pyabc:distance_function: scale_type not recognized.")

            if val == 0:
                # in practise, this case should be rare (if only for numeric
                # reasons, so setting the weight to 1 should be safe)
                w[key] = 1
            else:
                w[key] = 1 / val

        # normalize weights to have mean 1. This has just the effect that the
        # epsilon will decrease more smoothly, but is not important otherwise.
        mean_weight = statistics.mean(list(w.values()))
        for key in w.keys():
            w[key] /= mean_weight

        self.w[t] = w

        # logging
        df_logger.debug("update distance weights = {}".format(w))


def median_absolute_deviation(data: List):
    """
    Calculate the sample `median absolute deviation (MAD)
    <https://en.wikipedia.org/wiki/Median_absolute_deviation/>`_, defined as
    median(abs(data - median(data)).

    Parameters
    ----------

    data: List
        List of data points.

    Returns
    -------

    mad
        The median absolute deviation of the data.

    """

    data_median = statistics.median(data)
    normalized_data = []
    for item in data:
        normalized_data.append(abs(item - data_median))
    mad = statistics.median(normalized_data)
    return mad


def standard_deviation(data: List):
    """
    Calculate the sample `standard deviation (SD)
    <https://en.wikipedia.org/wiki/Standard_deviation/>`_.

    Parameters
    ----------

    data: List
        List of data points.

    Returns
    -------

    sd
        The standard deviation of the data points.
    """

    sd = statistics.stdev(data)
    return sd
