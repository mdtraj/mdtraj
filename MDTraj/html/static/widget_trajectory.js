// require([
//     "widgets/js/widget",
//     "rmol",
//     ],
// function(WidgetManager, RMol) {
require([
    "widgets/js/widget",
    ],
function(WidgetManager) {
    var TrajectoryView = IPython.DOMWidgetView.extend({
        render : function() {
            var container = $("<div/>").css({
                            width: '302px',
                            height: '302px', border:'1px solid #ccc',
                            margin: 'auto'}).append(
                                $('<div id="inner">').css({
                                    position: 'absolute',
                                    width: '300px',
                                    height: '300px',
                                })
                        );
            this.setElement(container);
            
            console.log(container);
            this.setElement(container);
            console.log("Creating RMol");
            this.rmol = new RMol(container);
            this.rmol.enableMouse();
            this.update();            
        },
        
        update : function () {
            console.log('update, xyz[0]=', this.model.attributes._xyz[0]);

            this.rmol.setTopology(this.model.attributes._topology);        
            this.rmol.setXYZ(this.model.attributes._xyz);
            
            var representation = {
                color: 'chainbow',
                mainChain: 'thickRibbon',
                sideChains: 'line'
            };

            this.rmol.setRepresentation(representation);
            this.rmol.zoomInto(this.rmol.getAllAtoms());
            this.rmol.render();
            
            return TrajectoryView.__super__.update.apply(this);
        }
    });
    WidgetManager.register_widget_view('TrajectoryView', TrajectoryView);
});