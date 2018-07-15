from .odorant_list import BaseOdorantListView
from compounds.models import Odorant, UserNotes


class ChemFilterListView(BaseOdorantListView):
    paginate_by = 100

    def get_context_data(self, **kwargs):
        context = super(ChemFilterListView, self).get_context_data(**kwargs)
        context['chem_type'] = self.kwargs['chem_type']
        context['page_header'] = self.kwargs['chem_type'].replace('_', ' ')
        return context

    def get_queryset(self):
        unfiltered = super(ChemFilterListView, self).get_queryset
        filter_method = getattr(Odorant.objects, self.kwargs['chem_type'], unfiltered)
        return filter_method()


class UserChemFilterListView(ChemFilterListView):
    template_name = 'odorants/user_odorant_list.html'
    context_object_name = 'compound_list'

    def get_queryset(self):
        notes_qs = UserNotes.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Odorant.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = super(UserChemFilterListView, self).get_queryset().filter(id__in=cpd_id_list)
        return queryset
