from django.shortcuts import get_object_or_404
from django.views.generic import ListView

from compounds.models import Bioactive, BioactiveCore, Enzyme
from compounds.views.mixins import BioactiveContentMixin, BioactiveSearchFilterMixin


class OligosaccharideListView(BioactiveContentMixin, BioactiveSearchFilterMixin, ListView):
    template_name = 'bioactives/cores_match_list.html'
    context_object_name = 'bioactive_list'
    bioactive_vals = None
    biocore = None
    category = 2
    paginate_by = 32

    def get_queryset(self):
        self.biocore = get_object_or_404(BioactiveCore, name='Oligosaccharides')
        return self.biocore.bioactives.all()

    def get_context_data(self, **kwargs):
        context = super(OligosaccharideListView, self).get_context_data(**kwargs)
        context.update({
            'page_header': self.biocore.name,
            'enzymes': Enzyme.objects.filter(category=1),
        })
        if self.request.GET.getlist('selected_bioactives'):
            id_list = [int(a) for a in self.request.GET.getlist('selected_bioactives') if a.isnumeric()]
            self.bioactive_vals = Bioactive.objects.filter(
                id__in=id_list
            ).values(
                'chemical_name',
                'chemical_properties',
                'cid_number',
                'cid_number_2',
                'iupac_name'
            )
        if self.bioactive_vals:
            context.update({
                'cid_numbers': [{'number': b['cid_number_2'] or b['cid_number'],
                                 'name': b['chemical_name'][:23] + '...' if len(b['chemical_name']) > 25
                                 else b['chemical_name'][:23] or b['iupac_name'][:32]}
                                for b in self.bioactive_vals],
            })
        return context
