from django.views.generic.edit import CreateView
from django.http import JsonResponse
from django.core.exceptions import ObjectDoesNotExist
import cirpy
import pubchempy as pcp

from compounds.models import Odorant
from compounds.forms import OdorantCreateForm


class BioactiveCreateView(CreateView):
    model = Odorant
    form_class = OdorantCreateForm
    template_name = 'odorants/create_odorant.html'

    def form_valid(self, form):
        form.instance.created_by = self.request.user.profile
        return super().form_valid(form)


def process_cas(request):
    data = {}
    cas_no = request.GET.get('cas_number')
    try:
        obj = Odorant.objects.get(cas_number__exact=cas_no)
        data['object_exists'] = obj.get_absolute_url()
        data['object_exists_name'] = obj.iupac_name
        return JsonResponse(data)
    except ObjectDoesNotExist:
        pass
    try:
        smiles = cirpy.query(cas_no, 'smiles')[0].value
        cid_no = pcp.get_compounds(smiles, 'smiles')[0].cid
        if smiles and cid_no:
            data['iupac_name'] = cirpy.Molecule(smiles).iupac_name
            data['structure_url'] = 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no)
            data.update({'smiles': smiles, 'hidden_cid': cid_no})
    except IndexError:
        data['error'] = 'No compound found for this CAS number'
    return JsonResponse(data)
