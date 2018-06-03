from django.views import generic
from django.shortcuts import get_object_or_404

from compounds.models import Compound, OdorType
from compounds.forms import CompoundFilter


class BaseCompoundListView(generic.ListView):
    queryset = Compound.objects.all().order_by('-trade_name', 'iupac_name')
    paginate_by = 25

    def get_context_data(self, **kwargs):
        context = super(BaseCompoundListView, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        compound_filter = CompoundFilter(self.request.GET, queryset=self.object_list)
        context['compound_filter'] = compound_filter
        return context


# ToDo: Move the above compound_filter to a mixin so dont have to use in e.g. all compound list view


class CompoundListView(BaseCompoundListView):
    pass


class OdorTypeCompoundListView(BaseCompoundListView):
    paginate_by = 25
    template_name = 'compounds/odor_compound_list.html'

    def get_queryset(self):
        self.odor_type = get_object_or_404(OdorType, term=self.kwargs['odor'])
        return Compound.objects.filter(odor_categories=self.odor_type)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['odor_type'] = self.odor_type
        return context


# Filter by custom model manager (e.g. phenolic, --- do aliphatic, aromatic, heterocyclic etc. as in book

class PhenolListView(BaseCompoundListView):
    queryset = Compound.objects.all_phenols().order_by('-trade_name', 'iupac_name')
    paginate_by = 40
    template_name = 'compounds/phenol_list2.html'


class PhenolListView(BaseCompoundListView):
    queryset = Compound.objects.all_phenols().order_by('-trade_name', 'iupac_name')
    paginate_by = 40
    template_name = 'compounds/phenol_list2.html'





