from django.shortcuts import redirect, reverse
from django.utils.text import slugify
from django.views.generic import ListView

from compounds.models import Activity, Bioactive
from compounds.views.bioactive.bioactive_list import BaseBioactiveListView
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
    classification = None

    def get_context_data(self, **kwargs):
        context = super(BaseBioactiveActivityFilterListView, self).get_context_data(**kwargs)
        context['body_systems'] = [a[1] for a in Activity.classifications]
        context['drug_actions'] = Activity.objects.actions().order_by('name')
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
        return context


class BioactiveDrugActionListView(BaseBioactiveActivityFilterListView):
    action_id = None

    def dispatch(self, request, *args, **kwargs):
        actions = Activity.objects.filter(category=1).values_list('name', 'id')
        slug_map = dict([(slugify(a[0]), a[1]) for a in actions])
        self.action_id = slug_map[kwargs['action']]
        self.queryset = Bioactive.objects.filter(activity__action_id=self.action_id)
        if not self.queryset.exists():
            self.queryset = Bioactive.objects.filter(activity_id=self.action_id)
        return super(BioactiveDrugActionListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveDrugActionListView, self).get_context_data(**kwargs)
        context['current_action'] = self.kwargs['action']
        context['page_header'] = self.kwargs['action'].replace('-', ' ')
        context['mechanisms'] = Activity.objects.mechanisms().filter(action_id=self.action_id)
        return context


class BioactiveMechanismListView(BaseBioactiveActivityFilterListView):
    template_name = 'bioactives/bioactive_list.html'
    action_id = None

    def dispatch(self, request, *args, **kwargs):
        actions = Activity.objects.filter(category=1).values_list('name', 'id')
        action_slug_map = dict([(slugify(a[0]), a[1]) for a in actions])
        self.action_id = action_slug_map[kwargs['action']]
        mechanisms = Activity.objects.filter(action_id=self.action_id).values_list('name', 'id')
        slug_map = dict([(slugify(a[0]), a[1]) for a in mechanisms])
        self.queryset = Bioactive.objects.filter(activity_id=slug_map[kwargs['mechanism']])
        return super(BioactiveMechanismListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveMechanismListView, self).get_context_data(**kwargs)
        context.update({
            'current_action': self.kwargs['action'],
            'page_header': self.deslug_mechanism_title(self.kwargs['mechanism']),
            'mechanisms': Activity.objects.mechanisms().filter(action_id=self.action_id),
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
