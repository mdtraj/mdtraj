try:
    from .install import enable_notebook
    from .trajectory_widget import TrajectoryView
    from .linked_heatmap import TrajectoryHeatmap
except ImportError:
    pass
