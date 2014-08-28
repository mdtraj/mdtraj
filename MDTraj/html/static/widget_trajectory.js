require([
    "widgets/js/widget",
    "rmol",
    ],
function(WidgetManager, RMol) {
    var TrajectoryView = IPython.DOMWidgetView.extend({
        render : function() {
            var container = $("<div/>").css({
                            width: '502px',
                            height: '502px', border:'1px solid #ccc',
                            margin: 'auto'}).append(
                                $('<div id="inner">').css({
                                    position: 'absolute',
                                    width: '500px',
                                    height: '500px',
                                })
                        );
            this.setElement(container);
            
            console.log(container);
            this.setElement(container);
            console.log("Creating RMol");
            this.rmol = new RMol(container);
            window.model = this.model;
            this.rmol.setTopology(this.model.attributes._topology);

            representation = {
                color: 'chainbow',
                mainChain: 'thickRibbon',
                sideChains: 'line'
            };

            this.rmol.setRepresentation(representation);
            this.enableMouse();
            
            this.update();
            
        },
        
        update : function () {
            this.rmol.setXYZ(this.model.attributes._xyz);
            this.rmol.zoomInto(this.rmol.getAllAtoms());
            this.rmol.render();
        }
    });
    WidgetManager.register_widget_view('TrajectoryView', TrajectoryView);
});