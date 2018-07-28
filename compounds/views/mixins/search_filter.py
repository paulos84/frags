import re

from django.shortcuts import redirect, reverse

from compounds.forms import BioactiveSearchForm, OdorantSearchForm


class BioactiveSearchFilterMixin(object):

    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z0-9]')
        cas_no = request.GET.get('cas_number')
        if cas_no or request.GET.get('iupac_name'):
            return redirect(reverse(
                'odorant-name-filter',
                kwargs={'search_query': regex.sub('', cas_no) if cas_no
                        else regex.sub('', request.GET.get('iupac_name').lower())})
            )
        return super(BioactiveSearchFilterMixin, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(BioactiveSearchFilterMixin, self).get_context_data(**kwargs)
        context['compound_search'] = BioactiveSearchForm()
        return context


class OdorantSearchFilterMixin(object):

    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z0-9]')
        cas_no = request.GET.get('cas_number')
        if cas_no or request.GET.get('iupac_name'):
            return redirect(reverse(
                'odorant-name-filter',
                kwargs={'search_query': regex.sub('', cas_no) if cas_no
                        else regex.sub('', request.GET.get('iupac_name').lower())})
            )
        return super(OdorantSearchFilterMixin, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(OdorantSearchFilterMixin, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        return context
