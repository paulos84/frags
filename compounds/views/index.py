from django.views.generic import TemplateView
from django.shortcuts import render
from django.shortcuts import redirect, reverse

from compounds.forms import ChemNameSearchForm, ProteinSearchForm
from compounds.models import Activity, Bioactive, Odorant, OdorType


class IndexView(TemplateView):
    template_name = 'index.html'

    def get(self, request, *args, **kwargs):
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
        context.update({
            'mols_count': Bioactive.objects.count() + Odorant.objects.count(),
            'body_systems': Activity.classified_actions_mechs(),
            'odor_categories': OdorType.objects.all(),
            'compound_search': ChemNameSearchForm,
            'search_form': ProteinSearchForm,
        })
        return context


def handler404(request, exception):
    return render(request, 'registration/error_404.html', {})


def handler500(request, exception):
    return render(request, 'registration/error_500.html', {})
