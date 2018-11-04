from django.views.generic import ListView

from compounds.models import Activity, Bioactive
from compounds.views.mixins import BioactiveSearchFilterMixin


class BioactiveApprovalsListView(BioactiveSearchFilterMixin, ListView):
    paginate_by = 32
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'
    queryset = Bioactive.objects.filter(
        approval_date__isnull=False
    ).order_by(
        '-approval_date'
    )

    def get_context_data(self, **kwargs):
        context = super(BioactiveApprovalsListView, self).get_context_data(**kwargs)
        context.update({
            'approvals_list': True,
            'body_systems': [a[1] for a in Activity.classifications],
            'drug_actions': Activity.objects.actions().order_by('name'),
            'page_header': 'Recent FDA Approvals'
        })
        return context
