from django.views import generic
import pubchempy as pcp

from compounds.models import Compound


class CompoundDetailView(generic.DetailView):
    model = Compound

    def get_context_data(self, **kwargs):
        context = super(CompoundDetailView, self).get_context_data(**kwargs)
        odor_types = self.get_object().odor_categories.values_list('term')
        context['odor_types'] = ', '.join([a[0] for a in odor_types])
        cid_no = self.get_object().cid_number
        context['structure_url'] = self.get_object().structure_url
        try:
            context['synonyms'] = ', '.join(pcp.get_compounds(cid_no)[0].synonyms)
            # PARSE OUT EC... from synonyms and FEMA...
        except KeyError:
            context['synonyms'] = 'n/a'
        return context

    # override get_object  to add the select_related for ...  prefetch_related
