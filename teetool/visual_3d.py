# functions to visualise the information
# (trajectories / probability) in 3 dimensions

import numpy as np
from scipy.interpolate import griddata
import mayavi.mlab as mlab

import teetool as tt


class Visual_3d(object):
    """
    <description>
    """

    def __init__(self, thisWorld):
        """
        <description>
        """

        # start figure
        self._mfig = mlab.figure(size=(800,600))
        self._world = thisWorld

    def plotTrajectories(self, list_clusters):
        """
        <description>
        """

        colours = tt.helpers.getDistinctColours(len(list_clusters))

        for (i, icluster) in enumerate(list_clusters):
            this_cluster = self._world.getCluster(icluster)
            for (x, Y) in this_cluster["data"]:
                mlab.plot3d(Y[:, 0], Y[:, 1], Y[:, 2], color=colours[i],
                            tube_radius=None)

    def plotLogDifference(self, icluster1, icluster2):
        """
        plots difference
        """

        [xx, yy, zz] = self._world.getGrid(ndim=3,
                                           resolution=[50, 50, 50])

        ss = np.zeros_like(xx)

        # add
        this_cluster = self._world.getCluster(icluster1)
        if ("logp" in this_cluster):
            (Y, s) = this_cluster["logp"]
            s_min = np.min(s)
            # interpolate result
            ss1 = griddata(Y, s, (xx, yy, zz),
                           method='linear',
                           fill_value=s_min)

            ss += ss1

        # subtract
        this_cluster = self._world.getCluster(icluster2)
        if ("logp" in this_cluster):
            (Y, s) = this_cluster["logp"]
            s_min = np.min(s)
            # interpolate result
            ss1 = griddata(Y, s, (xx, yy, zz),
                           method='linear',
                           fill_value=s_min)

            ss -= ss1

        # normalise
        ss_norm = (ss - np.min(ss)) / (np.max(ss) - np.min(ss))

        # mayavi
        src = mlab.pipeline.scalar_field(xx, yy, zz, ss_norm)

        # plot a volume
        # mlab.pipeline.volume(src, vmin=pmin, vmax=pmax)
        # slice it
        mlab.pipeline.image_plane_widget(src,
                                         plane_orientation='z_axes',
                                         slice_index=10,
                                         )

    def plotTube(self, list_clusters, popacity=0.3):
        """
        plots log-probability
        """

        [xx, yy, zz] = self._world.getGrid(ndim=3,
                             resolution=[40, 40, 40])

        # ss = np.zeros_like(xx)
        nclusters = len(list_clusters)
        lcolours = tt.helpers.getDistinctColours(nclusters)

        for (i, icluster) in enumerate(list_clusters):
         this_cluster = self._world.getCluster(icluster)
         if ("tube" in this_cluster):
             (Y, s) = this_cluster["tube"]
             # interpolate result
             ss = griddata(Y, s, (xx, yy, zz),
                         method='linear',
                         fill_value=np.min(s))

             ss = np.squeeze(ss)  # TODO fix this at a previous step
             # mayavi
             src = mlab.pipeline.scalar_field(xx, yy, zz, ss)

             # plot an iso surface
             mlab.pipeline.iso_surface(src,
                                       contours=[0.5],
                                       opacity=popacity,
                                       color=lcolours[i])

    def plotLogIsoSurface(self, list_clusters, pcontours=[.1, .2], popacity=0.3):
        """
        plots log-probability
        """

        [xx, yy, zz] = self._world.getGrid(ndim=3,
                                resolution=[40, 40, 40])

        # ss = np.zeros_like(xx)
        nclusters = len(list_clusters)
        lcolours = tt.helpers.getDistinctColours(nclusters)

        for (i, icluster) in enumerate(list_clusters):
            this_cluster = self._world.getCluster(icluster)
            if ("logp" in this_cluster):
                (Y, s) = this_cluster["logp"]
                # interpolate result
                ss = griddata(Y, s, (xx, yy, zz),
                            method='linear',
                            fill_value=np.min(s))

                # normalise
                ss_norm = (ss - np.min(ss)) / (np.max(ss) - np.min(ss))

                # mayavi
                src = mlab.pipeline.scalar_field(xx, yy, zz, ss_norm)

                # plot an iso surface
                mlab.pipeline.iso_surface(src,
                                          contours=pcontours,
                                          opacity=popacity,
                                          color=lcolours[i])

    def plotLogProbability(self, list_clusters, pmin=0.0, pmax=1.0):
        """
        plots log-probability
        """

        [xx, yy, zz] = self._world.getGrid(ndim=3,
                                           resolution=[30, 30, 30])

        ss = np.zeros_like(xx)

        for icluster in list_clusters:
            this_cluster = self._world.getCluster(icluster)
            if ("logp" in this_cluster):
                (Y, s) = this_cluster["logp"]
                s_min = np.min(s)
                # interpolate result
                ss1 = griddata(Y, s, (xx, yy, zz),
                               method='linear',
                               fill_value=s_min)
                # sum
                ss += ss1

        # normalise
        ss_norm = (ss - np.min(ss)) / (np.max(ss) - np.min(ss))

        # mayavi
        src = mlab.pipeline.scalar_field(xx, yy, zz, ss_norm)

        # show peak areas
        #mlab.pipeline.iso_surface(src, contours=[0.5*s.ptp(), ], opacity=0.1)
        # plot a volume
        #mlab.pipeline.volume(src, vmin=pmin, vmax=pmax)
        # slice it
        mlab.pipeline.image_plane_widget(src,
                                         plane_orientation='z_axes',
                                         slice_index=10,
                                         opacity=0.5,
                                         vmin=pmin,
                                         vmax=pmax)

    def plotOutline(self):
        """
        adds an outline
        """

        plot_outline = self._world.getExpandedOutline()

        mlab.outline(extent=plot_outline)

    def _plotTitle(self):
        """
        adds a title
        """

        # add title
        world_name = self._world.getName()

        if not (world_name == None):
            mlab.title(world_name)

    def save(self, add=None):
        """
        saves as file
        """

        if (add==None):
            saveas = self._world.getName()
        else:
            saveas = "{0}_{1}".format(self._world.getName(), add)

        mlab.savefig("output/3d_{0}.png".format(saveas), figure=self._mfig)

    def show(self):
        """
        shows the image [waits for user input]
        """

        # show figure
        mlab.show()

    def close(self):
        """
        closes figure(s)
        """

        mlab.close(all=True)
