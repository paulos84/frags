from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.palettes import Spectral6
from bokeh.plotting import figure
from django.views.generic import TemplateView

from compounds.models import Odorant, Substructure, OdorType
from compounds.forms import ChemDataChoiceForm
from compounds.utils.general import chemical_properties_label_map
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class SubstructureListView(OdorantSearchFilterMixin, TemplateView):
    template_name = 'odorants/substructure_list.html'

    def get_context_data(self, **kwargs):
        context = super(SubstructureListView, self).get_context_data(**kwargs)
        context.update({
            'acyclic_terpenoids': Substructure.objects.acyclic_terpenoids(),
            'cyclic_terpenoids': Substructure.objects.cyclic_terpenoids(),
            'bicyclic_terpenoids': Substructure.objects.bicyclic_terpenoids(),
            'sesquiterpenoids': Substructure.objects.sesquiterpenoids(),
            'ionones_etc': Substructure.objects.ionones_etc(),
            'cyclic_ketones': Substructure.objects.cyclic_ketones(),
            'miscellaneous': Substructure.objects.miscellaneous(),
            'choice_form': ChemDataChoiceForm,
            'odor_types': OdorType.objects.values('term'),
        })
        property_choice = self.request.GET.get('property_choice')
        if property_choice:
            plot = self.make_plot(property_choice)
            script, div = components(plot, CDN)
            context['plot_script'] = script
            context['plot_div'] = div
        return context

    def make_plot(self, chem_property):
        averages = Substructure.compound_sets_averages(chem_property)
        plot_data = list(averages.keys()), list(averages.values())
        title = chemical_properties_label_map.get(chem_property, chem_property)
        source = ColumnDataSource(data=dict(substructures=plot_data[0], avg_vals=plot_data[1], color=Spectral6))
        max_val = max(plot_data[1])
        p = figure(x_range=plot_data[0], y_range=(0, max_val + max_val / 3), plot_height=250, title=title,
                   toolbar_location=None, tools="")
        p.vbar(x='substructures', top='avg_vals', width=0.9, color='color', source=source)
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p


class CompoundMatchSubstructureListView(SubstructureListView):
    template_name = 'odorants/odorant_substructures.html'
    compound = None

    def get_queryset(self):
        self.compound = Odorant.objects.get(id=self.kwargs['pk'])
        return Substructure.compound_matches(self.compound)

    def get_context_data(self, **kwargs):
        context = super(CompoundMatchSubstructureListView, self).get_context_data(**kwargs)
        context['compound'] = self.compound
        return context
