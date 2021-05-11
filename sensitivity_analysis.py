"""This file will generate boxplots of the variation of the trade-off weights"""
import pandas as pd
import matplotlib.pyplot as plt

# read data from csv format
data = pd.read_csv('sensitivity.csv', delimiter='\t').dropna()

# create plots
fig, ax = plt.subplots(figsize=(8, 6))
data.plot(kind='box', ax=ax)
ax.grid()
ax.set_ylabel('score')
plt.xticks(rotation=30)
fig.tight_layout()
plt.savefig(r"figures/sensitivity.pdf")