from django.views import generic

from compounds.models import Substructure, OdorType

from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.palettes import Spectral6
from bokeh.plotting import figure


class SubstructureListView(generic.ListView):
    queryset = Substructure.objects.all()
    template_name = 'compounds/substructure_list.html'

    def get_context_data(self, **kwargs):
        context = super(SubstructureListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        plot = self.make_plot()
        script, div = components(plot, CDN)
        context['the_script'] = script
        context['the_div'] = div
        return context

    def make_plot(self):
        fruits = ['Apples', 'Pears', 'Nectarines', 'Plums', 'Grapes', 'Strawberries']
        counts = [5, 3, 4, 2, 4, 6]
        source = ColumnDataSource(data=dict(fruits=fruits, counts=counts, color=Spectral6))
        p = figure(x_range=fruits, y_range=(0, 9), plot_height=250, title="Fruit Counts",
                   toolbar_location=None, tools="")
        p.vbar(x='fruits', top='counts', width=0.9, color='color', source=source)
        p.xgrid.grid_line_color = None
        return p

