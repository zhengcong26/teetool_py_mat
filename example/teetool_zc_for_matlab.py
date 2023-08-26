# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 10:27:03 2023

@author: cuilab
"""

# col = cluster_data[0][0]

# diff_col = np.diff(col)

import scipy.io as sio
import numpy as np
import teetool as tt

rotrj = sio.loadmat('2023-08-25_G_rotrj.mat')

rotrj = rotrj['rotrj']

test = rotrj[0, 0]

npoints = 10
x = np.linspace(-50, 50, num=npoints)

# In[2]:

for (k, test) in np.ndenumerate(rotrj):
    
    my_trajectories = []

    for (i, cell_input) in enumerate(test):
        Y = np.array(test[i, 0]).transpose()

        if np.isnan(Y).any():
            continue

        my_trajectories.append((x, Y))


# In[3]:

# create a world
    world = tt.World(name="test", ndim=2, resolution=[100, 100])

# add data
    world.addCluster(my_trajectories, "one")

# In[4]:

# visual
    visual = tt.visual_2d.Visual_2d(world)

# # plot 50 trajectories
# visual.plotTrajectories(ntraj=50)

# # obtain limits
# xlim = visual._ax.get_xlim()
# ylim = visual._ax.get_ylim()

# # set labels
# visual._ax.set_xlabel("x [m]")
# visual._ax.set_ylabel("y [m]")

# visual._ax.legend()


# In[5]:

# build the model
    settings = {"model_type": "resampling",
            "ngaus": 100}

    world.buildModel(settings)

# In[6]:

# visual
    visual = tt.visual_2d.Visual_2d(world)

    temp = visual.plotTube_data(sdwidth=2)

# # add confidence region(s)
# visual.plotTube(sdwidth=2)

# # set limits
# visual._ax.set_xlim(xlim)
# visual._ax.set_ylim(ylim)

# # set labels
# visual._ax.set_xlabel("x [m]")
# visual._ax.set_ylabel("y [m]")

    fname = f'rotrj_ci_{k}.mat'
    
# In[7]:
    sio.savemat(fname, {'z': temp[0], 'x': temp[1], 'y': temp[2]})
    
    
    
    
    
