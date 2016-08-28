# functions to visualise the information
# (trajectories / probability) in 2 dimensions

import numpy as np
import matplotlib.pyplot as plt

from teetool import helpers


class Visual_2d(object):
    """
    <description>
    """

    def __init__(self, thisWorld):
        """
        <description>
        """

        # start figure
        self._fig = plt.figure()
        self._ax = self._fig.gca()
        self._world = thisWorld

    def plotTrajectories(self, list_clusters):
        """
        <description>
        """

        colours = helpers.getDistinctColours(len(list_clusters))

        for (i, icluster) in enumerate(list_clusters):
            this_cluster = self._world.getCluster(icluster)
            for (x, Y) in this_cluster["data"]:
                self._ax.plot(Y[:, 0], Y[:, 1], color=colours[i])

    def plotLogProbability(self, list_clusters):
        """
        plots log-probability
        """

        [xx, yy] = self._world.getGrid()

        s = np.zeros_like(xx)

        for icluster in list_clusters:
            this_cluster = self._world.getCluster(icluster)
            if ("logp" in this_cluster):
                s += this_cluster["logp"]

        # normalise
        s = (s - np.min(s)) / (np.max(s) - np.min(s))

        # plot contours
        self._ax.contourf(xx, yy, s, 20)

    def plotOutline(self):
        """
        adds an outline
        """

        return True

    def show(self):
        """
        shows the image [waits for user input]
        """

        plt.show()
