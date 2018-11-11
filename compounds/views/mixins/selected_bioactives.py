from compounds.models import Bioactive


class SelectedBioactivesMixin:
    """
    Add values from Bioactive query to use in add context data upon get request
    """
    bioactive_vals = None

    def dispatch(self, request, *args, **kwargs):
        selected_bioactives = request.GET.getlist('selected_bioactives')
        if selected_bioactives:
            id_list = [int(a) for a in selected_bioactives if a.isnumeric()]
            self.bioactive_vals = Bioactive.objects.filter(
                id__in=id_list
            ).values(
                'chemical_name',
                'chemical_properties',
                'cid_number',
                'cid_number_2',
                'iupac_name'
            )
        return super(SelectedBioactivesMixin, self).dispatch(request, *args, **kwargs)
