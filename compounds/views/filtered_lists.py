from .compound_list import BaseCompoundListView
from compounds.models import Compound, CompoundNotes


class ChemFilterListView(BaseCompoundListView):
    queryset = Compound.objects.aliphatic_carbonyls()
    paginate_by = 100

    def get_context_data(self, **kwargs):
        context = super(ChemFilterListView, self).get_context_data(**kwargs)
        context['chem_type'] = self.kwargs['chem_type']
        context['page_header'] = self.kwargs['chem_type'].replace('_', ' ')
        return context

    def get_queryset(self):
        unfiltered = super(ChemFilterListView, self).get_queryset
        filter_method = getattr(Compound.objects, self.kwargs['chem_type'], unfiltered)
        return filter_method()


class UserChemFilterListView(ChemFilterListView):
    template_name = 'compounds/user_compound_list.html'
    context_object_name = 'compound_list'

    def get_queryset(self):
        notes_qs = CompoundNotes.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Compound.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = super(UserChemFilterListView, self).get_queryset().filter(id__in=cpd_id_list)
        return queryset
