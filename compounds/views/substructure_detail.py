from django.db.models import Q
from django.views.generic import ListView
from django.views.generic.detail import SingleObjectMixin

from compounds.models import Compound, UserNotes, OdorType, Substructure
from compounds.forms import CompoundSearchForm


class SubstructureDetail(SingleObjectMixin, ListView):
    paginate_by = 14
    template_name = "compounds/substructure_detail.html"

    def get(self, request, *args, **kwargs):
        self.object = self.get_object(queryset=Substructure.objects.all())
        return super(SubstructureDetail, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(SubstructureDetail, self).get_context_data(**kwargs)
        context['odor_types'] = OdorType.objects.values('term')
        context['compound_search'] = CompoundSearchForm()
        return context

    def get_queryset(self):
        qs = self.object.compound_set()
        cas_number = self.request.GET.get('cas_number')
        iupac_name = self.request.GET.get('iupac_name')
        if cas_number:
            qs = qs.filter(iupac_name__exact=cas_number)
        elif iupac_name:
            qs = qs.filter(Q(iupac_name__icontains=iupac_name) |
                           Q(trade_name__icontains=iupac_name) |
                           Q(chemical_properties__synonyms__icontains=iupac_name))
        return qs


class ChemFilterSubstructureDetail(SubstructureDetail):
    def get_queryset(self):
        unfiltered = super(ChemFilterSubstructureDetail, self).get_queryset
        filter_method = getattr(Compound.objects, self.kwargs['chem_type'], unfiltered)
        return Compound.substructure_matches(self.object.smiles, filter_method())

    def get_context_data(self, **kwargs):
        context = super(ChemFilterSubstructureDetail, self).get_context_data(**kwargs)
        context['chem_type'] = self.kwargs['chem_type'].replace('_', ' ')
        return context


class UserSubstructureDetail(SubstructureDetail):
    def get_queryset(self):
        notes_qs = UserNotes.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Compound.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Compound.objects.filter(id__in=cpd_id_list)
        return Compound.substructure_matches(self.object.smiles, queryset)
