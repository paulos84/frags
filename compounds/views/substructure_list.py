from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.palettes import Spectral6
from bokeh.plotting import figure
from django.views import generic

from compounds.models import Substructure, OdorType
from compounds.forms import ChemDataChoiceForm
from compounds.utils.general import chemical_properties_label_map


class SubstructureListView(generic.ListView):
    queryset = Substructure.objects.all()
    template_name = 'compounds/substructure_list.html'

    def get_context_data(self, **kwargs):
        context = super(SubstructureListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        context['choice_form'] = ChemDataChoiceForm
        plot = self.make_plot()
        script, div = components(plot, CDN)
        context['plot_script'] = script
        context['plot_div'] = div
        return context

    def make_plot(self, chem_property='mw'):
        averages = Substructure.compound_sets_averages('mw')
        plot_data = list(averages.keys()), list(averages.values())
        title = chemical_properties_label_map.get(chem_property, chem_property)
        source = ColumnDataSource(data=dict(substructures=plot_data[0], avg_vals=plot_data[1], color=Spectral6))
        p = figure(x_range=plot_data[0], y_range=(0, max(plot_data[1]) + 60), plot_height=250, title=title,
                   toolbar_location=None, tools="")
        p.vbar(x='substructures', top='avg_vals', width=0.9, color='color', source=source)
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p
