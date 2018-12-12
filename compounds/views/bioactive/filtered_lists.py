import json

from django.contrib import messages
from django.db.models import Q
from django.shortcuts import redirect, reverse
from django.utils.text import slugify
from django.views.generic import ListView

from compounds.models import Activity, Bioactive, Enzyme
from compounds.views.mixins import BioactiveSearchFilterMixin, SelectedBioactivesMixin


class BioactiveSearchFilterListView(BioactiveSearchFilterMixin, ListView):
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'
    queryset = Bioactive.objects.none()
    paginate_by = 32

    def dispatch(self, request, *args, **kwargs):
        search_query = kwargs.pop('search_query', '')
        field = kwargs.pop('field', '')
        if field == 'inchikey':
            try:
                print(search_query)
                obj_id = Bioactive.objects.get(inchikey=search_query).id
                return redirect(reverse('bioactive-detail', kwargs={'pk': obj_id}))
            except Bioactive.DoesNotExist:
                messages.info(request, 'No compound matching InChiKey, please try adding it to the database')
                return redirect(reverse('bioactive-list', kwargs={'category': 'medicinal'}))
        if field == 'name':
            specific_matches = Bioactive.objects.filter(chemical_name__iexact=search_query)
            if specific_matches.exists():
                return redirect(reverse('bioactive-detail', kwargs={'pk': specific_matches.first().id}))
            chem_name_matches = Bioactive.objects.filter(
                Q(chemical_name__icontains=search_query) | Q(chemical_properties__synonyms__icontains=search_query)
                | Q(iupac_name__icontains=search_query))
            if chem_name_matches.count() == 1:
                return redirect(reverse('bioactive-detail', kwargs={'pk': chem_name_matches.first().id}))
            search_terms = search_query.split(' ')
            if chem_name_matches.exists():
                self.queryset = chem_name_matches
            elif len(search_terms) >= 2:
                messages.info(request, "No matches for: '{}'. Showing matches for: '{}'".format(
                    ' '.join(search_terms), search_terms[0]))
                self.queryset = Bioactive.objects.filter(chemical_name__icontains=search_terms[0])
            if not self.queryset.exists():
                messages.info(request, 'No compounds matched the search query')
                return redirect(reverse('bioactive-list', kwargs={'category': 'medicinal'}))
        else:
            self.queryset = Bioactive.objects.filter(iupac_name__icontains=search_query)
            if not self.queryset.exists():
                messages.info(request, 'No compounds matched the query')
                return redirect(reverse('bioactive-list', kwargs={'category': 'medicinal'}))
        return super(BioactiveSearchFilterListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveSearchFilterListView, self).get_context_data(**kwargs)
        context.update({
            'body_systems': Activity.classified_actions_mechs(),
            'drug_actions': Activity.objects.actions().order_by('name'),
        })
        return context


class BaseBioactiveActivityFilterListView(BioactiveSearchFilterMixin, ListView):
    model = Bioactive
    paginate_by = 32
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'
    classification = None

    def get_context_data(self, **kwargs):
        context = super(BaseBioactiveActivityFilterListView, self).get_context_data(**kwargs)
        context.update({
            'body_systems': Activity.classified_actions_mechs(),
            'drug_actions': Activity.objects.actions().order_by('name'),
        })
        return context


class BioactiveClassificationListView(BaseBioactiveActivityFilterListView):
    action = None

    def dispatch(self, request, *args, **kwargs):
        classification_map = {slugify(k): v for v, k in Activity.classifications}
        self.classification = classification_map[kwargs['classification']]
        self.queryset = Bioactive.objects.filter(activity__classification=self.classification)
        return super(BioactiveClassificationListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveClassificationListView, self).get_context_data(**kwargs)
        context['page_header'] = dict(Activity.classifications)[self.classification]
        context.update({
            'page_header': dict(Activity.classifications)[self.classification],
        })
        return context


class BioactiveDrugActionListView(SelectedBioactivesMixin, BaseBioactiveActivityFilterListView):
    action_id = None
    show_proteins = False

    def dispatch(self, request, *args, **kwargs):
        slug_map = Activity.slug_map()
        self.action_id = slug_map[kwargs['action']]
        self.queryset = Bioactive.objects.filter(activity__action_id=self.action_id)
        if not self.queryset.exists():
            self.queryset = Bioactive.objects.filter(activity_id=self.action_id)
        self.show_proteins = kwargs.pop('show_proteins', '')
        return super(BioactiveDrugActionListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveDrugActionListView, self).get_context_data(**kwargs)
        mechanisms = Activity.objects.mechanisms().filter(action_id=self.action_id)
        proteins = Enzyme.objects.filter(mechanism__in=mechanisms)
        context.update({
            'current_action':  self.kwargs['action'],
            'page_header': self.kwargs['action'].replace('-', ' '),
            'mechanisms': mechanisms,
            'd3_struct': True,
            'show_proteins': 'true' if self.show_proteins else 'false',
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


class BioactiveMechanismListView(SelectedBioactivesMixin, BaseBioactiveActivityFilterListView):
    template_name = 'bioactives/bioactive_list.html'
    action_id = None
    mechanism_id = None

    def dispatch(self, request, *args, **kwargs):
        actions = Activity.objects.filter(category=1).values_list('name', 'id')
        action_slug_map = dict([(slugify(a[0]), a[1]) for a in actions])
        self.action_id = action_slug_map[kwargs['action']]
        mechanisms = Activity.objects.filter(action_id=self.action_id).values_list('name', 'id')
        slug_map = dict([(slugify(a[0]), a[1]) for a in mechanisms])
        self.mechanism_id = slug_map[kwargs['mechanism']]
        self.queryset = Bioactive.objects.filter(activity_id=self.mechanism_id)
        return super(BioactiveMechanismListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveMechanismListView, self).get_context_data(**kwargs)
        proteins = Enzyme.objects.filter(mechanism_id=self.mechanism_id)
        context.update({
            'current_action': self.kwargs['action'],
            'page_header': self.deslug_mechanism_title(self.kwargs['mechanism']),
            'd3_struct': True,
            'mech_pk': self.mechanism_id,
            'proteins': proteins.exists(),
            'proteins_json': json.dumps(
                [{'citation': a['citation'], 'notes': a['notes'], 'pdb_number': a['pdb_number']}
                 for a in proteins.values('pdb_number', 'notes', 'citation')])
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

    @staticmethod
    def deslug_mechanism_title(page_header):
        page_header = page_header.strip('-').replace('-', ' ')
        deslug_map = {
            'ace': 'ACE',
            'adrenaline': 'β-adrenaline',
            'hmg coa': 'HMG-CoA',
            'cox': 'COX',
            'β2-adrenergic': 'β2-adrenergic',
            'α-adrenergic': 'α-adrenergic',
            'pde5': 'PDE5',
            'sglt 2': 'SGLT-2',
            '5 ht3': '5-HT3',
            'nk1': 'NK1',
        }
        for k in deslug_map.keys():
            if page_header.startswith(k):
                return page_header.replace(k, deslug_map[k])
        for a in ['β', 'μ', ]:
            if page_header.startswith(a):
                return page_header
        return page_header.capitalize()
