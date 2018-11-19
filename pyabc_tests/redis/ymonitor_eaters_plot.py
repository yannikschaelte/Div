import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


vals = []

with open("ymonitor_eaters_RES_" + str(sys.argv[1]) + ".log") as f:
    content = f.readlines()

content = [x.strip() for x in content]

# basic unit: kb

for x in content:
    if x.endswith('m'):
        val = int(float(x[:-1]) * 1e3)
    elif x.endswith('g'):
        val = int(float(x[:-1]) * 1e6)
    elif x.endswith('t'):
        val = int(float(x[:-1]) * 1e9)
    else:
        try:
            val = int(x)
        except Exception:
            continue
    # print(x, val)
    vals.append(val)

vals = np.array(vals, dtype=int)
max_val = max(vals)
times = np.linspace(1, len(vals), len(vals)) / 3600

plt.plot(times, np.log10(vals))
plt.xlabel("Time [h]")
plt.ylabel("log10(RES [MB])")
plt.title("Max: " + str(max_val))
plt.savefig("ymonitor_eaters_RES_" + str(sys.argv[1]) + ".png")
