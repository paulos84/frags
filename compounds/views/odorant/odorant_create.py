from django.views.generic.edit import CreateView
from django.http import JsonResponse
from django.core.exceptions import ObjectDoesNotExist
import cirpy
import pubchempy as pcp

from compounds.models import Odorant
from compounds.forms import OdorantCreateForm, OdorantSearchForm
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class OdorantCreateView(OdorantSearchFilterMixin, CreateView):
    model = Odorant
    form_class = OdorantCreateForm
    template_name = 'odorants/create_odorant.html'

    def form_valid(self, form):
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super(OdorantCreateView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        return context


def process_cas(request):
    cas_no = request.GET.get('cas_number')
    try:
        obj = Odorant.objects.get(cas_number__exact=cas_no)
        data = {
            'object_exists': obj.get_absolute_url(),
            'object_exists_name': obj.chemical_name,
        }
        return JsonResponse(data)
    except ObjectDoesNotExist:
        pass
    try:
        smiles = cirpy.query(cas_no, 'smiles')[0].value
        pcp_query = pcp.get_compounds(smiles, 'smiles')[0]
    except IndexError:
        return JsonResponse({
            'error': 'No compound found for this CAS number'
        })
    cid_no = pcp_query.cid
    if smiles and cid_no:
        data = {
            'chemical_name': Odorant.scrape_compound_name(cid_no),
            'iupac_name': pcp_query.iupac_name,
            'structure_url': 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no),
            'hidden_cid': cid_no,
            'smiles': smiles,
        }
        return JsonResponse(data)
