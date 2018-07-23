import re

from django.shortcuts import redirect, reverse

from compounds.forms import OdorantSearchForm


class SearchFilterMixin(object):

    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z0-9]')
        if request.GET.get('cas_number') or request.GET.get('iupac_name'):
            return redirect(reverse(
                'odorant-name-filter',
                kwargs={'search_query': regex.sub('', request.GET.get('cas_number')) if request.GET.get('cas_number')
                        else regex.sub('', request.GET.get('iupac_name').lower())})
            )
        return super(SearchFilterMixin, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(SearchFilterMixin, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        return context
