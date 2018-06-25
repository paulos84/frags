from django.views.generic import ListView
from django.views.generic.detail import SingleObjectMixin
from compounds.models import Compound, UserNotes, OdorType, Substructure


class SubstructureDetail(SingleObjectMixin, ListView):
    paginate_by = 14
    template_name = "compounds/substructure_detail.html"

    def get(self, request, *args, **kwargs):
        self.object = self.get_object(queryset=Substructure.objects.all())
        return super(SubstructureDetail, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(SubstructureDetail, self).get_context_data(**kwargs)
        context['substructure'] = self.object
        return context

    def get_queryset(self):
        return Compound.substructure_matches(self.object.smiles)


class ChemFilterSubstructureDetail(SubstructureDetail):
    def get_queryset(self):
        unfiltered = super(ChemFilterSubstructureDetail, self).get_queryset
        filter_method = getattr(Compound.objects, self.kwargs['chem_type'], unfiltered)
        return Compound.substructure_matches(self.object.smiles, filter_method())


class UserSubstructureDetail(SubstructureDetail):
    def get_queryset(self):
        notes_qs = UserNotes.objects.filter(user=self.request.user.profile).values('compound')
        if not notes_qs:
            return Compound.objects.none()
        cpd_id_list = [a['compound'] for a in notes_qs]
        queryset = Compound.objects.filter(id__in=cpd_id_list)
        return Compound.substructure_matches(self.object.smiles, queryset)
