import numpy as np
import scipy as sp


class YDict(dict):
	__getattr__ = dict.__getitem__
	__setattr__ = dict.__setitem__
	

def dhc(fun, x0, lb, ub, options=None):
	"""
	Optimization via Dynamic Hill Climbing. Based on 
	[De La Maza, Yuret. Dynamic Hill Climbing].

	Parameters
	----------
	fun : callable
		The objective function to be minimized.
			``fun(x) -> floag``
		where x is an 1-D array with shape(n,1).
	x0: ndarray, shape (n,)
		Initial guess. Array of real elements of size (n,),
		where 'n' is the number of independent variables.
	lb: ndarray, shape (n,)
		Lower bound.
	ub: ndarray, shape (n,)
		Upper bound.
	options: dict, optional
		A dictionary of solver options.
		
	Returns
	-------
	res : dict
		The optimization results represented as a dict. 
		Important keys are: ``x`` the solution array, 
		``fval`` the found function value.
	"""
	
	# interpret options
	options = dhc_options(options)
	
	# number of variables
	dim = len(x0)
	
	# running variables
	x = x0.copy
	fval = fun(x)
	
	return YDict(x=x,fval=fval)

	
def dhc_options(options_in=None) -> YDict:
	"""
	Create default options dictionary.
	
	Parameters
	----------
	options: dict, optional
		A dictionary with some options pre-set. Every key
		is checked for validity.
		
	Returns
	-------
	res: dict
		The optimization results.
	"""
	
	# initialize options with default values
	
	options = YDict()
	options.TolX                = 1e-8;
	options.TolFun              = 1e-8;
	options.MaxFunEvals         = np.inf;
	options.InitialStepSize     = 0.1;
	options.ExpandFactor        = 2.1;
	options.ContractFactor      = 0.47;
	options.StuckSearchFactor   = 4;
	options.Display             = 'off';
	
	# fill from input

	if options_in is None:
		options_in = YDict()
	
	for key in options_in:
		if not key in options:
			raise KeyError("Options field " + str(key) + " does not exist.")
		options[key] = options_in[key]
		
	return options