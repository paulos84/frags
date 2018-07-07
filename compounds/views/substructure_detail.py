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
        compound_filter = CompoundSearchForm(self.request.GET, queryset=Compound.objects.all())
        context['compound_filter'] = compound_filter
        return context

    def get_queryset(self):
        return Compound.substructure_matches(self.object.smiles) | Compound.iupac_name_matches(
            self.object.iupac_name_pattern)


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
