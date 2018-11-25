from django.contrib import messages
from django.shortcuts import redirect, reverse

from compounds.models import Odorant, UserOdorant
from compounds.views.odorant.odorant_list import BaseOdorantListView


class OdorantSearchFilterListView(BaseOdorantListView):

    def dispatch(self, request, *args, **kwargs):
        search_query = kwargs.pop('search_query', '').replace('_', ' ')
        field = kwargs.pop('field', '')
        if field == 'cas':
            try:
                obj_id = Odorant.objects.get(cas_number=search_query).id
                return redirect(reverse('odorant-detail', kwargs={'pk': obj_id}))
            except Odorant.DoesNotExist:
                messages.info(request, 'No compound matches found for the CAS number, try adding it to the database')
                return redirect(reverse('all-odorants'))
        if field == 'name':
            specific_matches = Odorant.objects.filter(chemical_name__iexact=search_query)
            if specific_matches.exists():
                return redirect(reverse('odorant-detail', kwargs={'pk': specific_matches.first().id}))
            chem_name_matches = Odorant.objects.filter(chemical_name__icontains=search_query)
            if chem_name_matches.count() == 1:
                return redirect(reverse('odorant-detail', kwargs={'pk': chem_name_matches.first().id}))
            search_words = search_query.split(' ')
            if chem_name_matches.exists():
                self.queryset = chem_name_matches
            elif len(search_words) >= 2:
                messages.info(request, "No matches for: '{}'. Showing matches for: '{}'".format(
                    ' '.join(search_words), search_words[0]))
                self.queryset = Odorant.objects.filter(chemical_name__icontains=search_words[0])
        if not self.queryset.exists():
            messages.info(request, 'No compound matching the search query, try adding it to the database')
            return redirect(reverse('all-odorants'))
        return super(OdorantSearchFilterListView, self).dispatch(request, *args, **kwargs)


class OdorantChemFilterListView(BaseOdorantListView):
    paginate_by = 100

    def get_context_data(self, **kwargs):
        context = super(OdorantChemFilterListView, self).get_context_data(**kwargs)
        context.update({
            'chem_type': self.kwargs['chem_type'],
            'page_header': self.kwargs['chem_type'].replace('_', ' '),
        })
        return context

    def get_queryset(self):
        if self.kwargs['chem_type'] == 'carbonyls':
            return Odorant.objects.filter(smiles__icontains='=O')
        elif self.kwargs['chem_type'] == 'alcohols':
            return Odorant.objects.filter(smiles__icontains='O').exclude(smiles__icontains='=O')
        elif self.kwargs['chem_type'] == 'heteroaromatics':
            return Odorant.objects.heteroaromatics()
        return self.queryset


class UserOdorantChemFilterListView(OdorantChemFilterListView):
    template_name = 'odorants/user_odorant_list.html'
    context_object_name = 'compound_list'

    def get_queryset(self):
        notes_qs = UserOdorant.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Odorant.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = super(UserOdorantChemFilterListView, self).get_queryset().filter(id__in=cpd_id_list)
        return queryset
