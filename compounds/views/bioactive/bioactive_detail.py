from django.urls import reverse
from django.views.generic import DetailView

from compounds.models import Bioactive


class BioactiveDetailView(DetailView):
    model = Bioactive
    template_name = 'bioactives/bioactive_detail.html'

    def get_context_data(self, **kwargs):
        context = super(BioactiveDetailView, self).get_context_data(**kwargs)
        chem_properties = self.get_object().chemical_properties
        key_map = {'mw': 'molecular_weight', 'hac': 'heavy atom count', 'hetac': 'heteroatom count',
                   'rbc': 'rotable_bond_count'}
        for key in key_map:
            if key in chem_properties:
                chem_properties[key_map[key]] = chem_properties.pop(key)
        context['chemical_properties'] = chem_properties
        return context

    def get_success_url(self):
        return reverse('bioactive-detail', kwargs={'pk': self.object.pk})
