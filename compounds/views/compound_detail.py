from django.views import generic
import pubchempy as pcp

from compounds.models.compound import Compound


class CompoundDetailView(generic.DetailView):
    model = Compound

    def get_context_data(self, **kwargs):
        context = super(CompoundDetailView, self).get_context_data(**kwargs)
        cid_no = self.get_object().cid_number
        context['structure_url'] = 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no)
        try:
            context['synonyms'] = ', '.join(pcp.get_compounds(cid_no)[0].synonyms)
        except KeyError:
            context['synonyms'] = 'n/a'
        return context