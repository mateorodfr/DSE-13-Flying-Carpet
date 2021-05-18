"""This file will do sensitivity analysis on the trade-off by dropping columns"""
import numpy as np
from itertools import chain, combinations
import matplotlib.pyplot as plt

# Get scoring arrays
concept1 = np.array([5, 2, 4, 2, 1, 2, 2])
concept2 = np.array([2, 1, 1, 1, 4, 1, 1])
concept3 = np.array([4, 4, 5, 3, 5, 3, 3])
concept4 = np.array([3, 5, 3, 4, 1, 4, 5])

trad_mat = np.vstack((concept1, concept2, concept3, concept4))

# Get weights array
weights = np.array([0.395, 0.23, 0.07, 0.039, 0.15, 0.053, 0.063])


def redistribute(weights, columnstodrop):
    """Takes the weights and redistributes the weight from the columns that will be dropped"""
    try:
        to_redistribute = np.sum(weights[columnstodrop])/(7 - len(columnstodrop))
        new_weights = weights + to_redistribute
        new_weights[columnstodrop] = 0
        return new_weights
    except RuntimeWarning:
        print(f"lenght of columns to drop too long: {columnstodrop}")


def powerset(iterable):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


# Summarise data
basescore = trad_mat @ weights.T
print(basescore)
scores = []
for columnstodrop in powerset(range(7)):
    coltodrop = list(columnstodrop)
    new_weights = redistribute(weights, coltodrop)
    # skip the one where we drop all columns
    if np.isclose(sum(new_weights), 1):
        score = trad_mat @ new_weights.T
        scores.append(list(score))

scores = np.array(scores)
# plotting

fig, ax = plt.subplots(figsize=(8, 6))

bp = ax.boxplot(scores, patch_artist=True)
for patch in bp['boxes']:
    patch.set(facecolor='pink')

ax.set_ylabel('Trade-off Score')
major_ticks = np.linspace(0.5, 5.5, 11)
minor_ticks = np.linspace(0.5, 5.5, 51)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)
ax.set_xticklabels(["Concept 1", "Concept 2", "Concept 3", "Concept 4"], rotation=30)

fig.tight_layout()
plt.show()
