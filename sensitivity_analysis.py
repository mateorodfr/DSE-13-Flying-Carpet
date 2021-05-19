"""This file will generate boxplots of the variation of the trade-off weights"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_boxplot(filename1, save=False):
    """Plots a boxplot based on trade-off table data"""
    # read data from csv format
    data = pd.read_csv(filename1, delimiter='\t').dropna()

    # create plots
    fig, ax = plt.subplots(figsize=(8, 6))

    bp = ax.boxplot(data, patch_artist=True)
    for patch in bp['boxes']:
        patch.set(facecolor='cyan')

    ax.set_ylabel('Trade-off Score')
    major_ticks = np.linspace(1.5, 4.5, 7)
    minor_ticks = np.linspace(1.5, 4.5, 31)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    ax.set_xticklabels(list(data.columns), rotation=30)
    fig.tight_layout()

    if save:
        plt.savefig(r"figures/sensitivity.pdf")
    else:
        plt.show()

if __name__ == "__main__":
    # plot_boxplot("tradeoff_data/sensitivity_combinations.csv")
    plot_boxplot("tradeoff_data/sensitivity_combinations.csv", save=True)
