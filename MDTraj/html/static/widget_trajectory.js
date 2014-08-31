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
            var container = $("<div/>").css({
                            width: (this.model.attributes.width + 2) + 'px',
                            height: (this.model.attributes.height + 2) + 'px',
                            border:'1px solid #ccc',
                            margin: 'auto'}).append(
                                $('<div id="inner">').css({
                                    position: 'absolute',
                                    width: this.model.attributes.width + 'px',
                                    height: this.model.attributes.height + 'px',
                                })
                        );
            // set the id attribute to something random
            var s4 = Math.floor((1 + Math.random()) * 0x10000).toString(16); 
            container.attr('id', s4);

            this.setElement(container);
            this.iview = new iview(container);
            this.update();            
        },
        
        update : function () {
            /* This could definitely be done more efficiently. Right now we're
            just recreating and redrawing everything. For the (presumably)
            common use case where you just want to update the positions to the
            next frame in a trajectory, there's no real need to redefine the
            topology and representation.
            */
            this.rmol.setTopology(this.model.attributes.topology);        
            this.rmol.setXYZ(this.model.attributes.coordinates);
            
            var representation = {
                color: this.model.attributes.color,
                mainChain: this.model.attributes.mainChain,
                sideChains: this.model.attributes.sideChains,
            };
            this.rmol.setRepresentation(representation);
            this.rmol.zoomInto(this.rmol.getAllAtoms());
            this.rmol.render();
            
            return TrajectoryView.__super__.update.apply(this);
        }
    });
    WidgetManager.register_widget_view('TrajectoryView', TrajectoryView);
});