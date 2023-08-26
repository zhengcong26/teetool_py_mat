#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import package
import teetool as tt


# In[2]:


# generate trajectory data

cluster_data_1 = tt.helpers.get_trajectories(ntype=0,
                                             ndim=2,
                                             ntraj=500,
                                             npoints=100,
                                             noise_std=0.0)

cluster_data_2 = tt.helpers.get_trajectories(ntype=1,
                                             ndim=2,
                                             ntraj=500,
                                             npoints=100,
                                             noise_std=0.0)


# In[3]:


# create a world
world = tt.World(name="toy", ndim=2, resolution=[100, 100])

# add data
world.addCluster(cluster_data_1, "one")
world.addCluster(cluster_data_2, "two")


# In[4]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

# plot 50 trajectories
visual.plotTrajectories(ntraj=50)

# obtain limits
xlim = visual._ax.get_xlim()
ylim = visual._ax.get_ylim()

# set labels
visual._ax.set_xlabel("x [m]")
visual._ax.set_ylabel("y [m]")

visual._ax.legend()

# save
# visual.save(add='1')

# show
# visual.show()


# In[5]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

# plot 50 trajectories
visual.plotTrajectories(ntraj=50)

# plot points along line
visual.plotTrajectoriesPoints(list_icluster=[0], x1=0.3, ntraj=50, marker='o', markersize=3)

# obtain limits
xlim = visual._ax.get_xlim()
ylim = visual._ax.get_ylim()

# set labels
visual._ax.set_xlabel("x [m]")
visual._ax.set_ylabel("y [m]")

visual._ax.legend()

# save
# visual.save(add='2')

# show
# visual.show()


# In[12]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

visual.plotTimeSeries(icluster=0, idim=1, ntraj=50, colour=(0.4, 0.1, 0.1))

# set labels
visual._ax.set_xlabel(r"$\tau$ [-]")
visual._ax.set_ylabel("y [m]")
visual._ax.set_ylim(ylim)

line_label = visual._ax.plot([0.3,0.3], ylim,'--r')

visual.plotLegend()

# save
# visual.save(add='3')

# show
# visual.show()


# In[13]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

visual.plotTimeSeries(icluster=0, idim=0, ntraj=50, colour=(0.4, 0.1, 0.1))

# set labels
visual._ax.set_xlabel(r"$\tau$ [-]")
visual._ax.set_ylabel("x [m]")
visual._ax.set_ylim(xlim)

line_label = visual._ax.plot([0.3,0.3], xlim,'--r')

visual.plotLegend()

# save
# visual.save(add='4')

# show
# visual.show()


# In[8]:


# build the model
settings = {"model_type":"resampling",
            "ngaus":100}

world.buildModel(settings)


# In[9]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

# add confidence region(s)
visual.plotTube(sdwidth=2)

# plot points along line
visual.plotTrajectoriesPoints(list_icluster=[0], x1=0.3, ntraj=50, marker='o', markersize=3)

# set limits
visual._ax.set_xlim(xlim)
visual._ax.set_ylim(ylim)

# set labels
visual._ax.set_xlabel("x [m]")
visual._ax.set_ylabel("y [m]")

# save
# visual.save(add='5')

# show
# visual.show()


# In[10]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

# add confidence region(s)
visual.plotTube(sdwidth=2)

# add complexity map
cax = visual.plotComplexityMap(list_icluster=None, complexity=1, resolution=[100, 100], cmap1='hot_r')
# add colourbar
cbar = visual.plotColourBar(cax, ticks=[0, 1], label='Probability', shrink=0.7)
cbar.ax.set_yticklabels(['Low', 'High'])

# set limits
visual._ax.set_xlim(xlim)
visual._ax.set_ylim(ylim)

# set labels
visual._ax.set_xlabel("x [m]")
visual._ax.set_ylabel("y [m]")

# save
# visual.save(add='6')

# show
# visual.show()


# In[11]:


# get_ipython().run_line_magic('', 'matplotlib inline')

# visual
visual = tt.visual_2d.Visual_2d(world, dpi=300)

# add confidence region(s)
visual.plotTube(sdwidth=2)

# add complexity map
cax = visual.plotComplexityMap(list_icluster=None, complexity=2, resolution=[100, 100], cmap1='hot_r')
# add colourbar
cbar = visual.plotColourBar(cax, ticks=[0, 1], label='Probability', shrink=0.7)
cbar.ax.set_yticklabels(['Low', 'High'])

# set limits
visual._ax.set_xlim(xlim)
visual._ax.set_ylim(ylim)

# set labels
visual._ax.set_xlabel("x [m]")
visual._ax.set_ylabel("y [m]")

# save
# visual.save(add='7')

# show
# visual.show()

