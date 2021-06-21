"""Script to make fancy plots from the budget breakdown"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def weight_distribution(masses, massfractions, mass_err, colors, tick_names):
    """Function to plot one budget breakdown"""
    x_arr = np.arange(len(masses))
    start_arr = np.cumsum(massfractions) - massfractions

    fig, (ax1, ax2) = plt.subplots(2, figsize=(7, 5), gridspec_kw={'height_ratios': [9, 1]})

    ax1.bar(x_arr, masses, yerr=mass_err, color=colors, edgecolor="black", capsize=5.0, tick_label=tick_names, alpha=0.7)
    ax1.set_ylabel("Mass of subsystem [kg]")

    for i, (width, color) in enumerate(zip(massfractions, colors)):
        ax2.barh([0], width, left=start_arr[i], color=color, alpha=0.7, height=0.1)

    ax2.set_xlabel("Mass fractions of subsystem [%]")
    ax2.set_yticks([])

    major_ticks = np.linspace(0, 600, 7)
    minor_ticks = np.linspace(0, 600, 31)

    ax1.set_yticks(major_ticks)
    ax1.set_yticks(minor_ticks, minor=True)

    ax1.grid(which='minor', alpha=0.2)
    ax1.grid(which='major', alpha=0.5)

    ax2.set_xlim(0, 100)
    ax2.set_xticks(np.linspace(0, 100, 21), minor=True)
    fig.tight_layout()
    return fig, ax1, ax2


def weight_distribution_comparison(results, err_dict, colors, tick_names):
    """Function to plot the comparison of two budget breakdowns"""
    fig, (ax1, ax2) = plt.subplots(2, figsize=(7, 5), gridspec_kw={'height_ratios': [9, 1]})
    labels = list(results.keys())
    for label, masses in results.items():
        x_arr = np.arange(len(masses))
        masses = np.array(masses)
        massfractions = masses / total_mass * 100
        start_arr = np.cumsum(massfractions) - massfractions

        if label == "Final":
            ax1.bar(x_arr, masses, yerr=err_dict[label], color=colors, edgecolor="red",
                    capsize=5.0, tick_label=tick_names, alpha=0.7, align="edge", width=0.4, lw=2.5, ls="solid")
        else:
            ax1.bar(x_arr, masses, yerr=err_dict[label], color=colors, edgecolor="black",
                    capsize=5.0, tick_label=tick_names, alpha=0.7, align="edge", width=-0.4, ls="solid", lw=2.5)

    ax1.set_ylabel("Mass of subsystem [kg]")

    for i, (width, color) in enumerate(zip(massfractions, colors)):
        ax2.barh([0], width, left=start_arr[i], color=color, alpha=0.8, height=0.1, edgecolor="white", lw=2)

    ax2.set_xlabel("Mass fractions of subsystems (final version) [%]")
    ax2.set_yticks([])

    major_ticks = np.linspace(0, 800, 9)
    minor_ticks = np.linspace(0, 800, 41)

    ax1.set_yticks(major_ticks)
    ax1.set_yticks(minor_ticks, minor=True)

    ax1.grid(which='minor', alpha=0.2)
    ax1.grid(which='major', alpha=0.5)
    ax1.set_axisbelow(True)

    ax2.set_xlim(0, 100)
    ax2.set_xticks(np.linspace(0, 100, 21), minor=True)
    fig.tight_layout()

    legend_elements = [Line2D([0], [0], color='black', lw=2.5, label=labels[0]),
                       Line2D([0], [0], color='red', lw=2.5, label=labels[1])]

    ax1.legend(handles=legend_elements, loc='upper right')

    return fig, ax1, ax2


def horizontal_barplot(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('jet')(
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(8, 2))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.4,
                label=colname, color=color, alpha=0.8)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            ax.text(x, y, str(int(c)), ha='center', va='center',
                    color=text_color)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')

    return fig, ax


if __name__ == "__main__":
    category_names = ['Payload', 'Batteries', 'Propulsion', 'Structure', 'Controller']
    results = {
        'Baseline': [600, 517, 160, 333.4, 56.6],
        'Final': [600, 353, 350.5, 219.345, 18]
    }

    std_dev = {
        "Baseline": [15, 50, 25, 10, 5],        # This line is in percent !!!
        "Final": [20, 50, 22, 40, 5.1]
    }
    key_i = "Baseline"
    std_dev[key_i] = [m*p/100 for m, p in zip(results[key_i], std_dev[key_i])]
    print("standard deviations for baseline review: ", sum(std_dev[key_i]))
    mass_estimation = np.array([600, 290 + 63, 270.5 + 80, 219.345, 18])
    total_mass = np.sum(mass_estimation)
    print("Total mass at Final review: ", total_mass)
    mass_percentages = mass_estimation / total_mass * 100
    mass_std = np.array([20, 50, 22, 40, 5.1])
    print("Total variance in mass Final review: ", np.sum(mass_std))

    category_colors = plt.get_cmap('jet')(np.linspace(0.15, 0.85, len(mass_estimation)))

    # horizontal_barplot(results, category_names)
    # weight_distribution(mass_estimation, mass_percentages, mass_std, category_colors, category_names)
    weight_distribution_comparison(results, std_dev, category_colors, category_names)
    # plt.savefig("figures/budget_breakdown_masses.pdf")
    plt.show()
