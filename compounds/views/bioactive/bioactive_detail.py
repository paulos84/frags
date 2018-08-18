from django.urls import reverse
from django.views.generic import DetailView

from compounds.models import Bioactive, BioactiveCore, UserBioactive
from compounds.views.mixins import BioactiveSearchFilterMixin


class BioactiveDetailView(BioactiveSearchFilterMixin, DetailView):
    model = Bioactive
    template_name = 'bioactives/bioactive_detail.html'

    def get_context_data(self, **kwargs):
        compound = self.get_object()
        context = super(BioactiveDetailView, self).get_context_data(**kwargs)
        chem_properties = compound.chemical_properties
        key_map = {'mw': 'molecular weight', 'hac': 'heavy atom count', 'hetac': 'heteroatom count',
                   'rbc': 'rotable bond count', 'bond_stereo_count': 'stereogenic bond count',
                   'h_bond_donor_count': 'H-bond donor count', 'h_bond_acceptor_count': 'H-bond acceptor count',
                   'atom_stereo_count': 'stereogenic atom count'}
        for key in key_map:
            if key in chem_properties:
                chem_properties[key_map[key]] = chem_properties.pop(key)
        chem_properties.pop('synonyms', None)
        if self.request.user.is_authenticated:
            try:
                user_compound = UserBioactive.objects.get(
                    compound=compound,
                    user=self.request.user.profile
                )
                context['user_data'] = user_compound.chemical_data
            except UserBioactive.DoesNotExist:
                pass
        context.update({
            'chemical_properties': chem_properties,
            'substructures': BioactiveCore.compound_matches(compound),
            'cid_string': compound.cid_number_2 or compound.cid_number,
        })
        return context

    def get_success_url(self):
        return reverse('bioactive-detail', kwargs={'pk': self.object.pk})
