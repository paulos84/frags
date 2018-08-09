from django.views.generic import ListView
from django.shortcuts import redirect, reverse
from django.utils.text import slugify

from compounds.models import Activity, Bioactive
from compounds.views.bioactive.bioactive_list import BaseBioactiveListView, BioactiveListView
from compounds.views.mixins import BioactiveSearchFilterMixin


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


class BaseBioactiveActivityFilterListView(BioactiveSearchFilterMixin, ListView):
    model = Bioactive
    paginate_by = 32
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'

    def get_context_data(self, **kwargs):
        context = super(BaseBioactiveActivityFilterListView, self).get_context_data(**kwargs)
        context['page_header'] = dict(Activity.classifications)[self.classification]
        context['body_systems'] = [a[1] for a in Activity.classifications]
        context['drug_actions'] = Activity.objects.actions()
        context['mechanisms'] = Activity.objects.mechanisms()
        return context


class BioactiveClassificationListView(BaseBioactiveActivityFilterListView):
    classification = None

    def dispatch(self, request, *args, **kwargs):
        classification_map = {slugify(k): v for v, k in Activity.classifications}
        self.classification = classification_map[kwargs['classification']]
        self.queryset = Bioactive.objects.filter(activity__classification=self.classification)
        return super(BioactiveClassificationListView, self).dispatch(request, *args, **kwargs)

class BioactiveDrugActionListView(BaseBioactiveActivityFilterListView):
    classification = None

    def dispatch(self, request, *args, **kwargs):
        classification_map = {slugify(k): v for v, k in Activity.classifications}
        self.classification = classification_map[kwargs['classification']]
        self.queryset = Bioactive.objects.filter(activity__classification=self.classification)
        return super(BioactiveClassificationListView, self).dispatch(request, *args, **kwargs)

class BioactiveMechanismListView(BaseBioactiveActivityFilterListView):
    classification = None

    def dispatch(self, request, *args, **kwargs):
        classification_map = {slugify(k): v for v, k in Activity.classifications}
        self.classification = classification_map[kwargs['classification']]
        self.queryset = Bioactive.objects.filter(activity__classification=self.classification)
        return super(BioactiveClassificationListView, self).dispatch(request, *args, **kwargs)

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
