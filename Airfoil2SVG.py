import os
import matplotlib.pyplot as plt
import numpy as np
import pprint

def main():

    theta = np.radians(-45)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))

    datFile = np.genfromtxt(r"C:\Users\marvd\Documents\GitHub\DSE-13-Flying-Carpet\airfoils\all_airfoils\coord_seligFmt\s8037.dat", skip_header=1, dtype = float)

    top = np.asarray([R.dot(v) for v in datFile[0:43]])
    bottom = np.asarray([R.dot(v) for v in datFile[43:]])

    pp = pprint.PrettyPrinter(indent=4)

    fig, ax = plt.subplots(1, figsize=(11, 11), dpi = 135)
    ax.plot(top[:, 0], top[:, 1], c="black", lw=3)
    ax.plot(bottom[:, 0], bottom[:, 1], c="black", lw=3)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.51, 1.51)
    ax.set_ylim(-0.91, 0.91)

    pos = np.asarray([0.25, 0.00])
    ax.arrow(R.dot(pos)[0], R.dot(pos)[1], 0, .75, length_includes_head = True, width = 0.005, ec = "Black", fc="Black", label = "T")
    ax.arrow(R.dot(pos)[0], R.dot(pos)[1], R.dot(pos)[0], .75, length_includes_head = True, width = 0.005, ec = "Black", fc="Black", label = "L")
    ax.arrow(R.dot(pos)[0], R.dot(pos)[1], 0.75, 0., length_includes_head = True, width = 0.005, ec = "Black", fc="Black", label = "Q")
    ax.arrow(R.dot(pos)[0], R.dot(pos)[1], 0.75, R.dot(pos)[1], length_includes_head = True, width = 0.005, ec = "Black", fc="Black", label = "D")
    plt.show()

    # fig.savefig('demo.svg', transparent=True)


if __name__ == "__main__":
    main()