import base64
import matplotlib.pyplot as pp
from six.moves import StringIO
from IPython.utils.traitlets import Any
from IPython.html.widgets import DOMWidget, CallbackDispatcher
from IPython.utils.traitlets import Unicode, CUnicode, Bool, Bytes
from matplotlib.backends.backend_agg import FigureCanvasAgg

__all__ = ['MPLFigureButton', 'ImageButton']


class ImageButton(DOMWidget):
    disabled = Bool(False, help="Enable or disable user changes.", sync=True)
    _view_name = Unicode('ImageButtonView', sync=True)

    format = Unicode('png', sync=True)
    width = CUnicode(sync=True)
    height = CUnicode(sync=True)
    _b64value = Unicode(sync=True)
    
    value = Bytes()
    def _value_changed(self, name, old, new):
        self._b64value = base64.b64encode(new)

    def __init__(self, **kwargs):
        super(ImageButton, self).__init__(**kwargs)
        self._click_handlers = CallbackDispatcher()
        self.on_msg(self._handle_button_msg)

    def on_click(self, callback, remove=False):
        self._click_handlers.register_callback(callback, remove=remove)
    
    def _handle_button_msg(self, _, content):
        if content.get('event', '') == 'click':
            self._click_handlers(self, content) 


class MPLFigureButton(ImageButton):
    fig = Any()
    
    def __init__(self, **kwargs):
        super(MPLFigureButton, self).__init__(**kwargs)
        if self.fig is None:
            self.fig = pp.gcf()
    
    def _fig_changed(self, name, old, new):
        self._sync_b64value(new)
        
    def _sync_b64value(self, figure):
        buf = StringIO()
        FigureCanvasAgg(figure).print_png(buf)
        self._b64value = base64.b64encode(buf.getvalue())
