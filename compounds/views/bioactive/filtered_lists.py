from django.shortcuts import redirect, reverse

from compounds.models import Activity, Bioactive
from compounds.views.bioactive.bioactive_list import BaseBioactiveListView, BioactiveListView


class BioactiveSearchFilterListView(BaseBioactiveListView):

    def dispatch(self, request, *args, **kwargs):
        search_query = kwargs.pop('search_query', '')
        field = kwargs.pop('field', '')
        if field == 'inchikey':
            try:
                obj_id = Bioactive.objects.get(inchikey=search_query).id
                return redirect(reverse('bioactive-detail', kwargs={'pk': obj_id}))
            except Bioactive.DoesNotExist:
                pass
        if field == 'name':
            specific_matches = Bioactive.objects.filter(chemical_name__iexact=search_query)
            if specific_matches.exists():
                return redirect(reverse('bioactive-detail', kwargs={'pk': specific_matches.first().id}))
            chem_name_matches = Bioactive.objects.filter(chemical_name__icontains=search_query)
            if chem_name_matches.count() == 1:
                return redirect(reverse('bioactive-detail', kwargs={'pk': chem_name_matches.first().id}))
            elif chem_name_matches.exists():
                self.queryset = chem_name_matches
        else:
            self.queryset = Bioactive.objects.filter(iupac_name__icontains=search_query)
        return super(BioactiveSearchFilterListView, self).dispatch(request, *args, **kwargs)


class BioactiveClassificationListView(BioactiveListView):
    classification = None

    def dispatch(self, request, *args, **kwargs):
        self.classification = Activity.objects.filterkwargs['classification']
        return super(BaseBioactiveListView, self).dispatch(request, *args, **kwargs)


#
# class BioactiveChemFilterListView(BaseBioactiveListView):
#     paginate_by = 100
#
#     def get_context_data(self, **kwargs):
#         context = super(BioactiveChemFilterListView, self).get_context_data(**kwargs)
#         context['chem_type'] = self.kwargs['chem_type']
#         context['page_header'] = self.kwargs['chem_type'].replace('_', ' ')
#         return context
#
#     def get_queryset(self):
#         unfiltered = super(BioactiveChemFilterListView, self).get_queryset
#         filter_method = getattr(Bioactive.objects, self.kwargs['chem_type'], unfiltered)
#         return filter_method()
#
#
# class UserBioactiveChemFilterListView(BioactiveChemFilterListView):
#     template_name = 'odorants/user_odorant_list.html'
#     context_object_name = 'compound_list'
#
#     def get_queryset(self):
#         notes_qs = UserBioactive.objects.filter(user=self.request.user.profile).values('compound')
#         if not notes_qs:
#             return Bioactive.objects.none()
#         cpd_id_list = [a['compound'] for a in notes_qs]
#         queryset = super(UserBioactiveChemFilterListView, self).get_queryset().filter(id__in=cpd_id_list)
#         return queryset
