require(["widgets/js/widget", "widgets/js/manager"], function(widget, manager){
    var ImageButtonView = widget.DOMWidgetView.extend({
        render : function(){
            // Called when view is rendered.
            this.setElement($("<img />"));
            this.update(); // Set defaults.
        },

        update : function(){
            // Update the contents of this view
            //
            // Called when the model is changed.  The model may have been 
            // changed by another view or by a state update from the back-end.
            var image_src = 'data:image/' + this.model.get('format') + ';base64,' + this.model.get('_b64value');
            console.log(this.model);
            this.$el.attr('src', image_src);

            var width = this.model.get('width');
            if (width !== undefined && width.length > 0) {
                this.$el.attr('width', width);
            } else {
                this.$el.removeAttr('width');
            }
            
            var height = this.model.get('height');
            if (height !== undefined && height.length > 0) {
                this.$el.attr('height', height);
            } else {
                this.$el.removeAttr('height');
            }
            return ImageButtonView.__super__.update.apply(this);
        },

        events: {
            // Dictionary of events and their handlers.
            'click': '_handle_click',
        },
        
        _handle_click: function(ev) {
            // Handles when the button is clicked.
            console.log(this.$el.offset());
            var top = this.$el.offset().top;
            var left = this.$el.offset().left;
            var xAspect = this.$el.width() / this.$el[0].naturalWidth;
            var yAspect = this.$el.height() / this.$el[0].naturalHeight;

            
            var x = (ev.pageX - left) / xAspect;
            var y = (this.$el.height() - (ev.pageY - top)) / yAspect;
            this.send({event: 'click', 'mouseX': x, 'mouseY': y});
        },
    });

    // Register the DatePickerView with the widget manager.
    manager.WidgetManager.register_widget_view('ImageButtonView', ImageButtonView);
});
