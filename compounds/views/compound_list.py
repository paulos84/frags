from django.views import generic
from django.shortcuts import get_object_or_404

from compounds.models import Compound, OdorType
from compounds.forms import CompoundFilter


class BaseCompoundListView(generic.ListView):
    queryset = Compound.objects.all().order_by('-trade_name', 'iupac_name')
    template_name = 'compounds/compound_list.html'
    paginate_by = 25

    def get_context_data(self, **kwargs):
        context = super(BaseCompoundListView, self).get_context_data(**kwargs)
        # context['odor_types'] = '  '.join([a.term for a in )
        compound_filter = CompoundFilter(self.request.GET, queryset=self.object_list)
        context['compound_filter'] = compound_filter
        return context


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

# Use a listviewset for the bolow?

class AliphaticCarbonylsListView(BaseCompoundListView):
    queryset = Compound.objects.aliphatic_carbonyls()


class AliphaticAlcoholsListView(BaseCompoundListView):
    queryset = Compound.objects.aliphatic_alcohols()


class AromaticCarbonylsListView(BaseCompoundListView):
    queryset = Compound.objects.aromatic_carbonyls()


class AromaticAlcoholsListView(BaseCompoundListView):
    queryset = Compound.objects.aromatic_alcohols()


class HeteroaromaticsListView(BaseCompoundListView):
    queryset = Compound.objects.heteroaromatics()