import re

from django.shortcuts import redirect, reverse

from compounds.forms import BioactiveSearchForm, OdorantSearchForm


class BioactiveSearchFilterMixin:
    """
    Enables the compound search form functionality by providing a method to handle GET requests as well as form itself
    """
    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z0-9]')
        chem_name = request.GET.get('chemical_name')
        iupac = request.GET.get('iupac_name', '').lower()
        inchikey = request.GET.get('inchikey', '')
        if any([inchikey, chem_name, iupac]):
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
        context['compound_search'] = BioactiveSearchForm()
        return context


class OdorantSearchFilterMixin:
    """
    Enables the compound search form functionality by providing a method to handle GET requests as well as form itself
    """
    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z0-9]')
        cas_no = request.GET.get('cas_number')
        chem_name = request.GET.get('chemical_name')
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
