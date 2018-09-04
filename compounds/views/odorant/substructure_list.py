from math import pi

from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from django.views.generic import TemplateView, ListView

from compounds.models import Odorant, Substructure, OdorType
from compounds.forms import ChemDataChoiceForm
from compounds.utils.chem_data import chemical_properties_label_map
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class SubstructureListView(OdorantSearchFilterMixin, TemplateView):
    template_name = 'odorants/substructure_list.html'

    def get_context_data(self, **kwargs):
        context = super(SubstructureListView, self).get_context_data(**kwargs)
        context.update({
            'substructure_sets': [
                {'subset': Substructure.objects.acyclic_terpenoids(), 'label': 'Acyclic Terpenoids'},
                {'subset': Substructure.objects.cyclic_terpenoids(), 'label': 'Cyclic Terpenoids'},
                {'subset': Substructure.objects.bicyclic_terpenoids(), 'label': 'Bicyclic Terpenoids'},
                {'subset': Substructure.objects.sesquiterpenoids(), 'label': 'Sespuiterpenoids'},
                {'subset': Substructure.objects.cycloaliphatic_ketones(), 'label': 'Damascones and Ionones'},
                {'subset': Substructure.objects.miscellaneous(), 'label': 'Miscellaneous'},
            ],
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

    @staticmethod
    def make_plot(chem_property):
        averages = Substructure.compound_sets_averages(chem_property)
        plot_data = list(averages.keys()), list(averages.values())
        title = chemical_properties_label_map.get(chem_property, chem_property)
        colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd',
                  '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
                  '#17becf', '#9edae5', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c', '#fdae6b',
                  '#fdd0a2', '#31a354', '#74c476', '#a1d99b', '#c7e9c0', '#756bb1', '#9e9ac8',	'#bcbddc', '#dadaeb',
                  '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
        source = ColumnDataSource(data=dict(substructures=plot_data[0], avg_vals=plot_data[1],
                                            color=colors[:len(plot_data[0])]))
        max_val = max(plot_data[1])
        p = figure(x_range=plot_data[0], y_range=(0, max_val + max_val / 3), plot_height=350, title=title,
                   toolbar_location=None, tools="")
        p.vbar(x='substructures', top='avg_vals', width=0.9, color='color', source=source)
        p.xaxis.major_label_orientation = pi / 4
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p


class CompoundMatchSubstructureListView(OdorantSearchFilterMixin, ListView):
    template_name = 'odorants/odorant_substructures.html'
    compound = None
    context_object_name = 'substructure_list'

    def dispatch(self, request, *args, **kwargs):
        self.compound = Odorant.objects.get(id=self.kwargs['pk'])
        return super(CompoundMatchSubstructureListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        return Substructure.compound_matches(self.compound)

    def get_context_data(self, **kwargs):
        context = super(CompoundMatchSubstructureListView, self).get_context_data(**kwargs)
        context['compound'] = self.compound
        return context
