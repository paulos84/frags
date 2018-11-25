from django.contrib import messages
from django.core.exceptions import ObjectDoesNotExist
from django.http import JsonResponse
from django.views.generic.edit import CreateView
import cirpy
import pubchempy as pcp

from compounds.models import Odorant
from compounds.forms import OdorantCreateForm, OdorantSearchForm
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class OdorantCreateView(OdorantSearchFilterMixin, CreateView):
    model = Odorant
    form_class = OdorantCreateForm
    template_name = 'odorants/create_odorant.html'

    def form_invalid(self, form):
        if 'iupac_name' in form.errors:
            messages.error(self.request, 'No valid data could be resolved for this compound')
        elif 'smiles' in form.errors:
            messages.error(self.request, 'SMILES string indicates this is not an odorant')
        return super().form_invalid(form)

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
            'object_exists_name': str(obj),
        }
        return JsonResponse(data)
    except ObjectDoesNotExist:
        pass
    try:
        smiles = cirpy.query(cas_no, 'smiles')[0].value
        pcp_query = pcp.get_compounds(smiles, 'smiles')[0]
        cid_no = pcp_query.cid
    except IndexError:
        return JsonResponse({
            'error': 'No compound found for this CAS number'
        })
    if smiles and cid_no:
        data = {
            'chemical_name': Odorant.scrape_compound_name(cid_no),
            'iupac_name': pcp_query.iupac_name,
            'structure_url': 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no),
            'hidden_cid': cid_no,
            'smiles': smiles,
        }
        return JsonResponse(data)
