"""Script to make fancy plots from the budget breakdown"""

import numpy as np
import matplotlib.pyplot as plt

# TODO refine data, right now it's just placeholders
# m_total = 1667
# # Payload, Batteries, structure, Propulsion, Controller
# subsys_mass = np.array([600, 617, 333.4, 160, np.nan])
# subsys_mass[-1] = m_total - np.sum(subsys_mass)
#
# subsys_perc = subsys_mass / m_total * 100
# subsys_err = np.array([2, 20, 10, 15, 5]) / 100
#
# x_arr = np.arange(len(subsys_perc))
#
# fig, axes = plt.subplots(2, figsize=(6, 8))
#
# ax[0].bar(x_arr, subsys_mass, yerr=subsys_err * subsys_mass)
#
# for i, width in enumerate(subsys_perc):
#     ax[1].barh([0], width, left=)
#
# fig.tight_layout()
# plt.show()


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


category_names = ['Payload', 'Batteries',
                  'Structure', 'Propulsion', 'Controller']
results = {
    'Baseline': [600, 517, 333.4, 160, 56.6],
    'Final': [600, 517, 333.4, 160, 56.6]
}

horizontal_barplot(results, category_names)
plt.show()
