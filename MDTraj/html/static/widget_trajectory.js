/*
This script creates and registers the TrajectoryView widget on the
browser side. Basically, it's a small div, with the RMol viewer (WebGL-based
protein visualization) hooked in. Changes to the class on the python side
propagate here and modify `this.model.attributes`, and re-call `update`.
*/

require([
    "widgets/js/widget",
    "iview",
    ],
function(WidgetManager, iview) {
    var TrajectoryView = IPython.DOMWidgetView.extend({
        render : function() {
            console.log('TrajectoryView.render');
            var $el = $("<canvas/>").css({
                            width: this.model.attributes.width + 'px',
                            height: this.model.attributes.height + 'px',
                        });
            this.setElement($el);
            this.iv = new iview($el);
            this.update();
        },

        update : function () {
            /* This could definitely be done more efficiently. Right now we're
            just recreating and redrawing everything. For the (presumably)
            common use case where you just want to update the positions to the
            next frame in a trajectory, there's no real need to redefine the
            topology and representation.
            */

            console.log('TrajectoryView.update');
            window.iv = this.iv;
            window.model = this.model;
            this.iv.loadTopology(this.model.attributes._topology);
            this.iv.loadCoordinates(this.model.attributes._coordinates);
            this.iv.loadSecondaryStructure(this.model.attributes._secondaryStructure);
            this.iv.zoomInto();
            return TrajectoryView.__super__.update.apply(this);
        }
    });
    WidgetManager.register_widget_view('TrajectoryView', TrajectoryView);
});