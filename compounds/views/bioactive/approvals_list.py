import json

from django.views.generic import ListView

from compounds.models import Activity, Bioactive, Enzyme
from compounds.views.mixins import BioactiveSearchFilterMixin, SelectedBioactivesMixin


class BioactiveApprovalsListView(SelectedBioactivesMixin, BioactiveSearchFilterMixin, ListView):
    paginate_by = 32
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'
    queryset = Bioactive.objects.filter(approval_date__isnull=False).order_by('-approval_date')

    def get_context_data(self, **kwargs):
        context = super(BioactiveApprovalsListView, self).get_context_data(**kwargs)
        id_list = []
        for a in self.queryset:
            if a.activity and a.activity.enzymes.exists() and a.activity.category == 2:
                for enz in a.activity.enzymes.all():
                    id_list.append(enz.id)
        proteins = Enzyme.objects.filter(pk__in=id_list)
        context.update({
            'approvals_list': True,
            'body_systems': Activity.classified_actions_mechs(),
            'drug_actions': Activity.objects.actions().order_by('name'),
            'page_header': 'Recent FDA Approvals',
            'd3_struct': True,
            'proteins': proteins.exists(),
            'proteins_json': json.dumps(
            [{'citation': a['citation'], 'notes': a['notes'], 'pdb_number': a['pdb_number']}
             for a in proteins.values('pdb_number', 'notes', 'citation')]),
        })
        if self.bioactive_vals:
            context.update({
                'data_display': 'true',
                'cid_numbers': [{'number': b['cid_number_2'] or b['cid_number'],
                                 'name': b['chemical_name'][:23] + '...' if len(b['chemical_name']) > 25
                                 else b['chemical_name'][:23] or b['iupac_name'][:32]}
                                for b in self.bioactive_vals],
            })
        return context
