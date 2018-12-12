import re

from django.shortcuts import redirect, reverse

from compounds.forms import BioactiveSearchForm, OdorantSearchForm, ProteinSearchForm

regex = re.compile('[^-a-z\sA-Z0-9_]')


class BioactiveSearchFilterMixin:
    """
    Enables the compound search form functionality by providing a method to handle GET requests as well as form itself
    """
    def get(self, request, *args, **kwargs):
        print(request.GET)
        protein_term = request.GET.get('protein_term')
        if protein_term:
            return redirect(reverse(
                'proteins', kwargs={'search_query': protein_term}
            ))
        chem_name = request.GET.get('chemical_name')
        iupac = request.GET.get('iupac_name', '').lower()
        inchikey = request.GET.get('inchikey', '')
        if any([inchikey, chem_name, iupac]):
            print(regex.sub('', inchikey))
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
        return super(BioactiveSearchFilterMixin, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveSearchFilterMixin, self).get_context_data(**kwargs)
        context.update({
            'compound_search': BioactiveSearchForm(),
            'protein_search': ProteinSearchForm(),
        })
        return context


class OdorantSearchFilterMixin:
    """
    Enables the compound search form functionality by providing a method to handle GET requests as well as form itself
    """
    def get(self, request, *args, **kwargs):
        cas_no = request.GET.get('cas_number', '').strip()
        chem_name = request.GET.get('chemical_name', '').replace(' ', '_').strip()
        iupac = request.GET.get('iupac_name', '').lower()
        if any([cas_no, chem_name, iupac]):
            params = (iupac, 'iupac')
            if cas_no:
                params = cas_no, 'cas'
            elif chem_name:
                params = chem_name, 'name'
            return redirect(reverse(
                'odorant-name-filter',
                kwargs={
                    'search_query': regex.sub('', params[0]),
                    'field': params[1],
                }
            ))
        return super(OdorantSearchFilterMixin, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(OdorantSearchFilterMixin, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        return context
