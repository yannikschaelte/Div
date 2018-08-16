import amici
import amici.plotting
import pesto
import importlib
import os
import sys
import numpy as np


model_name = 'model_conversion_reaction'
model_output_dir = model_name

# load amici module
sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()
model.setTimepoints(amici.DoubleVector(np.linspace(0, 10, 11)))
model.setParameterScale(amici.AMICI_SCALING_LOG10)
model.setParameters(amici.DoubleVector([-0.3, -0.7]))
solver = model.getSolver()
solver.setSensitivityMethod(amici.AMICI_SENSI_FSA)
solver.setSensitivityOrder(amici.AMICI_SENSI_ORDER_FIRST)
rdata = amici.runAmiciSimulation(model, solver, None)
edata = amici.ExpData(rdata['ptr'].get(), 1, 0)
#edata.setObservedDataStdDev(1.0)
# run from amici
#print(edata)
#print(amici.getExpDataFieldAsNumPyArray(edata, 'sigmay'))

rdata2 = amici.runAmiciSimulation(model, solver, edata)
for key in rdata2:
    print(key, ": ", rdata2[key])

# pesto stuff
objective = pesto.objective.AmiciObjective(model, solver, [edata], 1)
print(objective.get_fval([0, 0]))
