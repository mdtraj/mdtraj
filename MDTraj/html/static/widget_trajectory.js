require([
    "widgets/js/widget",
    "rmol",
    ],
function(WidgetManager, RMol) {
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
            this.setElement(container);
            this.rmol = new RMol(container);
            this.rmol.enableMouse();
            this.update();            
        },
        
        update : function () {
            this.rmol.setTopology(this.model.attributes.topology);        
            this.rmol.setXYZ(this.model.attributes.coordinates);
            
            var representation = {
                color: this.model.attributes.color,
                mainChain: this.model.attributes.mainChain,
                sideChains: this.model.attributes.sideChains,
            };
            console.log('representation', representation);

            this.rmol.setRepresentation(representation);
            this.rmol.zoomInto(this.rmol.getAllAtoms());
            this.rmol.render();
            
            return TrajectoryView.__super__.update.apply(this);
        }
    });
    WidgetManager.register_widget_view('TrajectoryView', TrajectoryView);
});