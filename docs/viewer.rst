IPython Notebook Viewer
=======================

.. This file needs to be updated for use with nglview. In the mean time, it has
   been removed from the toctree


Example
-------
This is an example of the trajectory viewer with PDBID 1F39. Running inside the IPython notebook as a widget,
many other options can be controled such as the color scheme, full screen, different representations, and more.
This `video <https://www.youtube.com/watch?v=Lwy2Hdsr518>`_ shows an example of it running inside the notebook.

.. raw:: html

   <div id="widgetview"></div>

   <script type="text/javascript" src="//cdnjs.cloudflare.com/ajax/libs/three.js/r68/three.min.js"></script>
   <script type="text/javascript" src="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.1/jquery-ui.min.js"></script>
   <script type="text/javascript" src="_static/trajectoryview-html/iview.js"></script>
   <script type="text/javascript" src="_static/trajectoryview-html/surface.min.js"></script>
   <script type="text/javascript" src="_static/trajectoryview-html/filesaver.js"></script>
   <script type="text/javascript" src="_static/trajectoryview-html/context.js"></script>

   <script>
   $.getJSON('_static/example-trajectoryview-data.json', function(data) {
       console.log(data);
       var HEIGHT = 300,
           WIDTH = 300,
           HEIGHT_PX = '300px',
           WIDTH_PX = '300px';       

       var canvas = $("<canvas/>").height(HEIGHT).width(WIDTH);
       var iv = new iview(canvas);
       var container = $('<div/>').css({width: HEIGHT_PX, height: WIDTH_PX})
           .resizable({
               aspectRatio: 1,
               resize: function(event, ui) {
                   iv.renderer.setSize(ui.size.width, ui.size.height);
                   },
	       stop : function(event, ui) {
	           iv.render()
	       },
           });
       container.append(canvas);

       $('#widgetview').append(container).css({
           'margin-left': 'auto',
           'margin-right': 'auto',
           'width': WIDTH
       });
       
       var options = {
           'camera': 'perspective',
           'background': 'white',
           'colorBy': 'spectrum',
           'primaryStructure': 'nothing',
           'secondaryStructure': 'cylinder & plate',
           'surface': 'nothing'
       };
        
       iv.loadTopology(data.topology);
       iv.loadCoordinates(data.frameData.coordinates);
       iv.loadAtomAttributes(data.frameData.secondaryStructure);
           
       iv.rebuildScene(options);
       iv.zoomInto(options);
       iv.render();


   });
   </script>

Usage
-----

To use the trajectory widget in the notebook, start your IPython session with ::

    from mdtraj.html import TrajectoryView, enable_notebook
    enable_notebook()

And then at the end of a cell, to render the viewer in the output space, use a ::

    TrajectoryView(your_mdtraj_trajectory)


Requirements
------------

    `IPython <http://ipython.org/>`_  >= 4.0
        To install IPython and the notebook with ``conda``,
        use ``conda install ipython-notebook ipywidgets``.

    Modern web browser with WebGL
        We've had the best luck with Chrome and modern versions
        of Firefox.
