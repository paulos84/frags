from math import pi

from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from django.shortcuts import get_object_or_404
from django.views.generic import ListView, TemplateView

from compounds.forms import ChemDataChoiceSubmitForm
from compounds.models import Activity, Bioactive, BioactiveCore
from compounds.utils.chem_data import chemical_properties_label_map, colors
from compounds.views.mixins.search_filter import BioactiveSearchFilterMixin


class BioactiveCoreListView(BioactiveSearchFilterMixin, TemplateView):
    template_name = 'bioactives/cores_list.html'

    def get_context_data(self, **kwargs):
        context = super(BioactiveCoreListView, self).get_context_data(**kwargs)
        context.update({
            'substructure_sets': [
                {'subset': BioactiveCore.objects.medicinal(), 'label': 'Medicinal'},
                {'subset': BioactiveCore.objects.food(), 'label': 'Nutraceutical'},
            ],
            'choice_form': ChemDataChoiceSubmitForm,
        })
        if self.request.GET.get('stats_data'):
            averages = BioactiveCore.compound_sets_stats()
            for chem_prop, val in averages.items():
                plot = self.make_plot(val, chem_prop, 'mean')
                if plot:
                    script, div = components(plot, CDN)
                    context['plot_script' + '_' + chem_prop] = script
                    context['plot_div' + '_' + chem_prop] = div
            stats_arrays = BioactiveCore.compound_sets_stats(std_dev=True)
            sds = {}
            for cp in chemical_properties_label_map.keys():
                sds[cp] = [(a[0], a[1][cp].std()) for a in stats_arrays]
            for chem_prop, val in sds.items():
                plot = self.make_plot(val, chem_prop, 'std dev in')
                if plot:
                    script, div = components(plot, CDN)
                    context['sd_plot_script' + '_' + chem_prop] = script
                    context['sd_plot_div' + '_' + chem_prop] = div
            context['data_display'] = 'true'
        return context

    @staticmethod
    def make_plot(stats_data, chem_property, title_description):
        plot_data = [a[0] for a in stats_data], [a[1] for a in stats_data]
        source = ColumnDataSource(data=dict(substructures=plot_data[0], avg_vals=plot_data[1],
                                            color=colors[:len(plot_data[0])]))
        max_val = max(plot_data[1])
        title = chemical_properties_label_map.get(
            chem_property, chem_property) + ' - ' + title_description + ' values for compounds of parent substructures'
        p = figure(x_range=plot_data[0], y_range=(0, max_val + max_val / 3), plot_height=350, title=title,
                   toolbar_location=None, tools="")
        p.vbar(x='substructures', top='avg_vals', width=0.9, color='color', source=source)
        p.xaxis.major_label_orientation = pi / 4
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p


class BioactiveCoreMatchList(BioactiveSearchFilterMixin, ListView):
    paginate_by = 28
    template_name = "bioactives/cores_match_list.html"
    model = BioactiveCore
    context_object_name = 'bioactive_list'
    bioactive_core = None
    bioactives = None
    names = []

    def dispatch(self, request, *args, **kwargs):
        self.bioactive_core = get_object_or_404(BioactiveCore, slug=kwargs['slug'])
        return super(BioactiveCoreMatchList, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        return self.bioactive_core.bioactive_set()

    def get_context_data(self, **kwargs):
        context = super(BioactiveCoreMatchList, self).get_context_data(**kwargs)
        context.update({
            'body_systems': [a[1] for a in Activity.classifications],
            'drug_actions': Activity.objects.actions().order_by('name'),
            'page_header': self.bioactive_core.name,
        })
        selected_bioactives = self.request.GET.getlist('selected_bioactives')
        if selected_bioactives:
            id_list = [int(a) for a in selected_bioactives if a.isnumeric()]
            bioactive_vals = Bioactive.objects.filter(
                id__in=id_list
            ).values(
                'chemical_name',
                'chemical_properties',
                'cid_number',
                'cid_number_2',
                'iupac_name'
            )
            self.bioactives = bioactive_vals
        if self.bioactives:
            context.update({
                'choice_form': ChemDataChoiceSubmitForm,
                'data_display': 'true',
                'cid_numbers': [{'number': b['cid_number_2'] or b['cid_number'],
                                 'name': b['chemical_name'][:23] + '...' if len(b['chemical_name']) > 25
                                 else b['chemical_name'][:23] or b['iupac_name'][:32]}
                                for b in self.bioactives],
            })
            properties = chemical_properties_label_map.keys()
            for cp in properties:
                plot = self.make_plot(cp)
                script, div = components(plot, CDN)
                context['plot_script' + '_' + cp] = script
                context['plot_div' + '_' + cp] = div
        return context

    def make_plot(self, chem_property):
        self.names = [b['chemical_name'][:23] + '...' if len(b['chemical_name']) > 25
                      else b['chemical_name'][:23] or b['iupac_name'][:32] for b in self.bioactives]
        raw_vals = [b['chemical_properties'][chem_property] for b in self.bioactives
                    if isinstance(b['chemical_properties']['synonyms'], str)]
        plot_data = self.names, [0 if a is None else a for a in raw_vals]
        title = chemical_properties_label_map.get(chem_property, chem_property)
        source = ColumnDataSource(data=dict(mechanisms=plot_data[0], avg_vals=plot_data[1],
                                            color=colors[:len(plot_data[0])]))
        max_val = max(plot_data[1])
        p = figure(x_range=plot_data[0], y_range=(0, max_val + max_val / 3), plot_height=350, title=title,
                   toolbar_location=None, tools="")
        p.vbar(x='mechanisms', top='avg_vals', width=0.9, color='color', source=source)
        p.xaxis.major_label_orientation = pi / 4
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p
