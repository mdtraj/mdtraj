require([
    "widgets/js/widget",
    "rmol",
    ],
function(WidgetManager, THREE) {
    console.log("Setting up TrajectoryView");
    console.log(THREE);
    console.log(THREE.TrackballControls);
    var TrajectoryView = IPython.DOMWidgetView.extend({
        render : function() {
            this.setElement($("<p>Hello !!@ World</p>"));
        }
    });
    WidgetManager.register_widget_view('TrajectoryView', TrajectoryView);
});