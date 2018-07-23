from django.views.generic import ListView
from django.db.models import Q

from compounds.forms import OdorantSearchForm
from compounds.views.mixins.search_filter import SearchFilterMixin


class BaseCompoundListView(SearchFilterMixin, ListView):
    paginate_by = 32

    def get_queryset(self):
        qs = super(BaseCompoundListView, self).get_queryset()
        cas_number = self.request.GET.get('cas_number')
        iupac_name = self.request.GET.get('iupac_name')
        if cas_number:
            qs = qs.filter(iupac_name__exact=cas_number)
        elif iupac_name:
            qs = qs.filter(Q(iupac_name__icontains=iupac_name) |
                           Q(trade_name__icontains=iupac_name) |
                           Q(chemical_properties__synonyms__icontains=iupac_name))
        return qs

    def get_context_data(self, **kwargs):
        context = super(BaseCompoundListView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        return context
