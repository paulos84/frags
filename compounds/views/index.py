import json

from django.http import JsonResponse
from django.views.generic import TemplateView
from django.shortcuts import render
from django.shortcuts import redirect, reverse

from compounds.forms import ChemNameSearchForm, ProteinSearchForm
from compounds.models import Activity, Bioactive, BioactiveCore, CompanyPipeline, Development, Odorant, OdorType


class IndexView(TemplateView):
    template_name = 'index.html'
    indications = Development.indication_pairs()
    classes = Development.classification_pairs()

    def get(self, request, *args, **kwargs):
        if request.is_ajax():
            return self.process_development(request)
        chem_name = request.GET.get('chemical_name')
        protein_term = request.GET.get('protein_term')
        if chem_name:
            return redirect(reverse(
                'bioactive-name-filter', kwargs={
                    'search_query': chem_name,
                    'field': 'name'}
            ))
        elif protein_term:
            return redirect(reverse(
                'proteins', kwargs={'search_query': protein_term}
            ))
        return super(IndexView, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(IndexView, self).get_context_data(**kwargs)
        company_vals = CompanyPipeline.objects.values('name', 'id')
        indications = set([(a[0], a[1]) for a in self.indications])
        classes = set([(a[0], a[1]) for a in self.classes])
        context.update({
            'mols_count': Bioactive.objects.count() + Odorant.objects.count(),
            'body_systems': Activity.classified_actions_mechs(),
            'biocores': BioactiveCore.objects.medicinal().values('name', 'cid_number', 'slug'),
            'odor_categories': OdorType.objects.all(),
            'compound_search': ChemNameSearchForm,
            'search_form': ProteinSearchForm,
            'companies': json.dumps(
                [{'id': a['id'], 'name': a['name']} for a in company_vals]),
            'indications': json.dumps(
                [{'id': a[0], 'name': a[1]} for a in indications]),
            'classifications': json.dumps(
                [{'id': a[0], 'name': a[1]} for a in classes if a[1]]),
            })
        return context

    def process_development(self, request):
        initial = request.GET.get('initial')
        action_id = request.GET.get('action_id', 0)
        class_choice = request.GET.get('class_choice')
        comp_choice_id = request.GET.get('company_id', 0)
        dev_cpds = []
        if initial:
            dev_cpds = Development.objects.all()
        elif comp_choice_id:
            dev_cpds = Development.objects.filter(company_id=int(comp_choice_id[0]))
        elif action_id:
            dev_ids = [a[2] for a in self.indications if a[0] == int(action_id)]
            dev_cpds = Development.objects.filter(id__in=dev_ids)
        elif class_choice:
            dev_cpds = Development.objects.filter(activity__classification=class_choice)
        dev_data = [
            {'phase': a.phase,
             'name': a.bioactive.chemical_name,
             'url': a.bioactive.get_absolute_url(),
             'structure_url': a.bioactive.structure_url,
             'activity': a.bioactive.activity.name if a.bioactive.activity else '',
             'activity_url': a.bioactive.activity_url,
             'company': a.company.name,
             'date': a.completion_date,
             'title': a.study_title,
             } for a in dev_cpds
            ]
        return JsonResponse({'dev_data': dev_data}, safe=False)


def handler404(request, exception):
    return render(request, 'registration/error_404.html', {})


def handler500(request, exception):
    return render(request, 'registration/error_500.html', {})
