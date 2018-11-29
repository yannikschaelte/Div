from model import *
import pyabc
import pyabc.visualization
import matplotlib.pyplot as plt


h = pyabc.History("sqlite:///db1.db")
h.id = 3
df, w = h.get_distribution(m=0, t=h.max_t)
ax = pyabc.visualization.plot_kde_1d(df, w, 'mean',
                                     xmin=limits['mean'][0], xmax=limits['mean'][1],
                                     numx=1000, refval=p_true)
plt.savefig("kde.png")
