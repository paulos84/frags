import json
import re

from django.contrib import messages
from django.db.models import Q
from django.shortcuts import redirect, reverse
from django.views.generic import TemplateView

from compounds.forms import BioactiveSearchForm, ProteinSearchForm
from compounds.models import Activity, Bioactive, Enzyme
from compounds.views.mixins import BioactiveSearchFilterMixin


class ActivityProteinListView(TemplateView):
    enz_match = False
    bioactivity_match = False
    search_term = None
    results_label = None
    template_name = 'bioactives/protein_list.html'

    def dispatch(self, request, *args, **kwargs):
        self.search_term = kwargs.pop('search_query', '')
        return super(ActivityProteinListView, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        protein_term = request.GET.get('search_term')
        if protein_term:
            return redirect(reverse(
                'proteins', kwargs={'search_query': protein_term}
            ))
        chem_name = request.GET.get('chemical_name')
        iupac = request.GET.get('iupac_name', '').lower()
        inchikey = request.GET.get('inchikey', '')
        if any([inchikey, chem_name, iupac]):
            regex = re.compile('[^-a-z\sA-Z0-9_]')
            params = (iupac, 'iupac')
            if inchikey:
                params = inchikey, 'inchikey'
            elif chem_name:
                params = chem_name, 'name'
            return redirect(reverse(
                'bioactive-name-filter',
                kwargs={
                    'search_query': regex.sub('', params[0]),
                    'field': params[1],
                }
            ))
        return super(ActivityProteinListView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(ActivityProteinListView, self).get_context_data(**kwargs)
        proteins = self.process_search_query(self.search_term)
        context.update({
            'body_systems': Activity.classified_actions_mechs(),
            'drug_actions': Activity.objects.actions().order_by('name'),
            'page_header': 'Protein search query: {}'.format(self.search_term),
            'search_form': ProteinSearchForm(),
            'compound_search': BioactiveSearchForm(),
            'results_label': ' {} drug targets'.format(self.results_label) if self.results_label else None,
            'proteins': proteins.exists(),
            'proteins_json': json.dumps(
                [{'citation': a['citation'], 'notes': a['notes'], 'pdb_number': a['pdb_number']}
                 for a in proteins.values('pdb_number', 'notes', 'citation')]),
        })
        return context

    def process_search_query(self, search_term):
        enzymes = Enzyme.objects.filter(notes__icontains=search_term)
        if enzymes.exists():
            return enzymes
        matches = Bioactive.objects.filter(
            Q(chemical_name__icontains=search_term) | Q(iupac_name__icontains=search_term))
        if not matches.exists():
            split_terms = search_term.split(' ')
            if len(split_terms) >= 2:
                matches = Bioactive.objects.filter(
                    Q(chemical_name__icontains=split_terms[0]) | Q(iupac_name__icontains=split_terms[0]))
        bioacts = [a['activity_id'] for a in matches.values('id', 'activity_id') if a['activity_id']]
        if bioacts:
            activity_id = max(bioacts)
            act_obj = Activity.objects.get(id=activity_id)
            action, self.results_label = (act_obj.action, act_obj.action.name) \
                if act_obj.category == 2 else (act_obj, act_obj.name)
            mechanism_ids = [a['id'] for a in action.mechanisms.values('id')]
            mech_enzs = Enzyme.objects.filter(mechanism__in=mechanism_ids)
            if mech_enzs.exists():
                return Enzyme.objects.filter(mechanism__in=mechanism_ids)
        messages.info(self.request, 'No proteins matched the query')
        return Enzyme.objects.exclude(notes__isnull=True)
