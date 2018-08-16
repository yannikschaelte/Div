import amici
import amici.plotting
import pesto

import importlib
import os
import sys
import numpy as np

sbml_file = 'model_conversion_reaction.xml'
model_name = 'model_conversion_reaction'
model_output_dir = model_name
sbml_importer = amici.SbmlImporter(sbml_file)
sbml_importer.sbml2amici(model_name, model_output_dir)
