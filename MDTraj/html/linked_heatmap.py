import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as pp
from IPython.html.widgets import ContainerWidget
from .trajectory_widget import TrajectoryView
from .imagebutton_widget import MPLFigureButton


def TrajectoryHeatmap(trajectory, x, y, fig=None, **kwargs):
    if fig is None:
        fig = pp.figure(figsize=(5,5))
        ax = fig.add_subplot(1,1,1)
        ax.hexbin(x, y)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    else:
        ax = fig.get_axes()[0]

    heatmap = MPLFigureButton(fig=fig, height='200px');
    viewer = TrajectoryView(trajectory, **kwargs)
    pointToDataCoordinates = ax.transData.inverted()
    kdtree = cKDTree(np.vstack((x, y)).T)

    def callback(b, content):
        x, y = pointToDataCoordinates.transform([content['mouseX'], content['mouseY']])
        _, index = kdtree.query(x=[x, y], k=1)
        print(x, y, index)
        viewer.frame = index

    heatmap.on_click(callback)
    return ContainerWidget(children=(heatmap, viewer))
