"""Script to make fancy plots from the budget breakdown"""

import numpy as np
import matplotlib.pyplot as plt


def weight_distribution(masses, massfractions, mass_err, colors, tick_names):
    x_arr = np.arange(len(masses))
    start_arr = np.cumsum(massfractions) - massfractions

    fig, (ax1, ax2) = plt.subplots(2, figsize=(7, 5), gridspec_kw={'height_ratios': [3, 1]})

    ax1.bar(x_arr, masses, yerr=mass_err, color=colors, tick_label=tick_names, alpha=0.5)
    ax1.set_ylabel("Mass of subsystem [kg]")

    for i, (width, color) in enumerate(zip(massfractions, colors)):
        ax2.barh([0], width, left=start_arr[i], color=color, alpha=0.5, height=0.6)

    ax2.set_xlabel("Mass fractions of subsystem [%]")
    ax2.set_yticks([])

    fig.tight_layout()
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
    category_colors = plt.get_cmap('RdYlGn')(
        np.linspace(0.15, 0.85, data.shape[1]))

    fig, ax = plt.subplots(figsize=(8, 2))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.4,
                label=colname, color=color)
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
    category_names = ['Payload', 'Batteries',
                      'Structure', 'Propulsion', 'Controller']
    results = {
        'Baseline': [600, 517, 333.4, 160, 56.6],
        'Final': [600, 517, 333.4, 160, 56.6]
    }

    # TODO: fix the masses to the correct subsystem values
    # TODO: Pick nicer colors
    # TODO: Add grid lines
    mass_estimation = np.array([600, 517, 333.4, 160, 56.6])
    mass_percentages = mass_estimation / np.sum(mass_estimation) * 100
    mass_std = np.array([20.3, 50, 22, 43, 5.1])

    category_colors = plt.get_cmap('Set2')(np.linspace(0, 1, len(mass_estimation)))

    horizontal_barplot(results, category_names)
    weight_distribution(mass_estimation, mass_percentages, mass_std, category_colors, category_names)

    plt.show()
