"""
<description>
"""

import pytest as pt
pt.importorskip("teetool.visual_2d")
import teetool as tt


def test_visual_2d():
    """
    produce figures
    """

    # from teetool import visual_2d

    # build world
    world_1 = tt.World(name="Example 3D", dimension=2)

    # extreme reduced resolution
    world_1.setResolution(xstep=3, ystep=3, zstep=3)

    # add trajectories
    for ntype in [0, 1]:
        correct_cluster_name = "toy {0}".format(ntype)
        correct_cluster_data = tt.helpers.get_trajectories(ntype, D=2, N=20)
        world_1.addCluster(correct_cluster_data, correct_cluster_name)

    # test grid
    [xx, yy] = world_1.getGrid()
    assert (xx.shape == yy.shape)

    # model all trajectories
    settings = {}
    settings["model_type"] = "resample"
    settings["mgaus"] = 10

    for i in [0, 1]:
        world_1.buildModel(i, settings)
        world_1.buildLogProbality(i)

    for i in [0, 1]:
        # visuals by mayavi
        visual = tt.visual_2d.Visual_2d(world_1)
        # visualise trajectories
        visual.plotTrajectories([i])
        # visualise intersection
        visual.plotLogProbability([i])
        # close
        visual.close()

    # visuals by mayavi
    visual = tt.visual_2d.Visual_2d(world_1)
    # visualise trajectories
    visual.plotTrajectories([0, 1])
    # visualise intersection
    visual.plotLogProbability([0, 1])
    # close
    visual.close()
