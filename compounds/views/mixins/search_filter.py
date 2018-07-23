import re

from django.shortcuts import redirect, reverse


class SearchFilterMixin(object):
    paginate_by = 32

    def get(self, request, *args, **kwargs):
        regex = re.compile('[^-a-z0-9]')
        if request.GET.get('cas_number') or request.GET.get('iupac_name'):
            return redirect(reverse(
                'odorant-name-filter',
                kwargs={'search_query': regex.sub('', request.GET.get('cas_number')) if request.GET.get('cas_number')
                        else regex.sub('', request.GET.get('iupac_name').lower())})
            )
        return super(SearchFilterMixin, self).get(request, *args, **kwargs)
