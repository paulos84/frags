from .compound_list import BaseCompoundListView
from compounds.models import Compound


class ChemFilterListView(BaseCompoundListView):
    queryset = Compound.objects.aliphatic_carbonyls()
    paginate_by = 100

    def get_context_data(self, **kwargs):
        context = super(ChemFilterListView, self).get_context_data(**kwargs)
        context['page_header'] = self.kwargs['chem_filter'].replace('_', ' ')
        return context

    def get_queryset(self):
        unfiltered = super(ChemFilterListView, self).get_queryset
        filter_method = getattr(Compound.objects, self.kwargs['chem_filter'], unfiltered)
        return filter_method()
