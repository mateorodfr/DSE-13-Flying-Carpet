"""This file will generate boxplots of the variation of the trade-off weights"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# read data from csv format
data = pd.read_csv('sensitivity.csv', delimiter='\t').dropna()

# create plots
fig, ax = plt.subplots(figsize=(8, 6))

bp = ax.boxplot(data, patch_artist=True)
for patch in bp['boxes']:
    patch.set(facecolor='cyan')

ax.set_ylabel('Trade-off Score')
major_ticks = np.linspace(2, 4, 5)
minor_ticks = np.linspace(2, 4, 25)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)
ax.set_xticklabels(list(data.columns), rotation=30)
fig.tight_layout()

# plt.show()
plt.savefig(r"figures/sensitivity.pdf")