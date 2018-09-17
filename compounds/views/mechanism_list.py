from math import pi

from bokeh.embed import components
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from bokeh.resources import CDN
from django.db.models import Count
from django.views.generic import ListView

from compounds.forms import ChemDataChoiceSubmitForm, ClassficationChoiceForm
from compounds.models import Activity
from compounds.utils.chem_data import chemical_properties_label_map, colors


class MechanismListView(ListView):
    model = Activity
    context_object_name = 'mechanism_list'
    template_name = 'bioactives/all_mechanisms.html'
    is_filtered = False

    def get_queryset(self):
        qs = Activity.objects.mechanisms().order_by(
            'action__name', 'name'
        ).annotate(
            item_count=Count('bioactives')
        )
        if self.request.GET.get('class_choice'):
            self.is_filtered = True
            return qs.filter(
               classification=self.request.GET['class_choice']
            )
        return qs

    def get_context_data(self, **kwargs):
        context = super(MechanismListView, self).get_context_data(**kwargs)
        context.update({
            'filter_form': ClassficationChoiceForm,
            'is_filtered': self.is_filtered,
        })
        selected_mechanisms = self.request.GET.getlist('selected_mechanisms')
        if selected_mechanisms:
            context.update({
                'choice_form': ChemDataChoiceSubmitForm,
                'data_display': 'true',
            })
            id_list = [int(a) for a in selected_mechanisms if a]
            if self.request.GET.get('mean_data'):
                averages = Activity.bioactives_data_stats(id_list)
                for chem_prop, val in averages.items():
                    plot = self.make_plot(val, chem_prop, 'mean')
                    if plot:
                        script, div = components(plot, CDN)
                        context['plot_script' + '_' + chem_prop] = script
                        context['plot_div' + '_' + chem_prop] = div
                stats_arrays = Activity.bioactives_data_stats(id_list, std_dev=True)
                sds = {}
                for cp in chemical_properties_label_map.keys():
                    sds[cp] = [(a[0], a[1][cp].std()) for a in stats_arrays]
                for chem_prop, val in sds.items():
                    plot = self.make_plot(val, chem_prop, 'std dev in')
                    if plot:
                        script, div = components(plot, CDN)
                        context['sd_plot_script' + '_' + chem_prop] = script
                        context['sd_plot_div' + '_' + chem_prop] = div
        return context

    @staticmethod
    def make_plot(stats_data, chem_property, title_description):
        plot_data = [a[0] for a in stats_data], [a[1] for a in stats_data]
        source = ColumnDataSource(data=dict(mechanisms=plot_data[0], avg_vals=plot_data[1],
                                            color=colors[:len(plot_data[0])]))
        max_val = max(plot_data[1])
        title = chemical_properties_label_map.get(
            chem_property, chem_property) + ' - ' + title_description + ' values for compounds with selected mechanisms'
        p = figure(x_range=plot_data[0], y_range=(0, max_val + max_val / 3), plot_height=350, title=title,
                   toolbar_location=None, tools="")
        p.vbar(x='mechanisms', top='avg_vals', width=0.9, color='color', source=source)
        p.xaxis.major_label_orientation = pi / 4
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p
