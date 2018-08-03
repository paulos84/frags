from math import pi

from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.palettes import Spectral6
from bokeh.plotting import figure
from django.shortcuts import get_object_or_404
from django.views.generic import ListView, TemplateView

from compounds.models import BioactiveCore
from compounds.views.mixins.search_filter import BioactiveSearchFilterMixin


class BioactiveCoreListView(BioactiveSearchFilterMixin, TemplateView):
    template_name = 'bioactives/cores_list.html'

    def get_context_data(self, **kwargs):
        context = super(BioactiveCoreListView, self).get_context_data(**kwargs)
        context.update({
            'substructure_sets': [
                {'subset': BioactiveCore.objects.medicinal(), 'label': 'Medicinal compounds'},
                {'subset': BioactiveCore.objects.phytochemical(), 'label': 'Phytochemicals'},
                {'subset': BioactiveCore.objects.miscellaneous(), 'label': 'Miscellaneous'},
            ],
        })
        # property_choice = self.request.GET.get('property_choice')
        # if property_choice:
        #     plot = self.make_plot(property_choice)
        #     script, div = components(plot, CDN)
        #     context['plot_script'] = script
        #     context['plot_div'] = div
        return context

    # @staticmethod
    # def make_plot(chem_property):
    #     averages = Bioactive.compound_sets_averages(chem_property)
    #     plot_data = list(averages.keys()), list(averages.values())
    #     title = chemical_properties_label_map.get(chem_property, chem_property)
    #     source = ColumnDataSource(data=dict(substructures=plot_data[0], avg_vals=plot_data[1], color=Spectral6))
    #     max_val = max(plot_data[1])
    #     p = figure(x_range=plot_data[0], y_range=(0, max_val + max_val / 3), plot_height=350, title=title,
    #                toolbar_location=None, tools="")
    #     p.vbar(x='substructures', top='avg_vals', width=0.9, color='color', source=source)
    #     p.xaxis.major_label_orientation = pi / 4
    #     p.xgrid.grid_line_color = None
    #     p.ygrid.grid_line_color = None
    #     return p


class BioactiveCoreMatchList(BioactiveSearchFilterMixin, ListView):
    paginate_by = 16
    template_name = "bioactives/bioactive_list.html"
    model = BioactiveCore
    context_object_name = 'bioactive_list'
    bioactive_core = None

    def dispatch(self, request, *args, **kwargs):
        self.bioactive_core = get_object_or_404(BioactiveCore, slug=kwargs['slug'])
        return super(BioactiveCoreMatchList, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        return self.bioactive_core.bioactive_set()

    def get_context_data(self, **kwargs):
        context = super(BioactiveCoreMatchList, self).get_context_data(**kwargs)
        context['page_header'] = self.bioactive_core.name
        return context
