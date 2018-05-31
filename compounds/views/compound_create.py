from django.views.generic.edit import CreateView
from django.http import JsonResponse
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Q
import cirpy
import pubchempy as pcp

from compounds.models import Compound
from compounds.forms.compound_forms import CompoundCreateForm


class CompoundCreateView(CreateView):
    model = Compound
    form_class = CompoundCreateForm
    template_name = 'compounds/create_compound.html'


def process_cas(request):
    cas_no = request.GET.get('cas_number')
    data = {}
    try:
        obj = Compound.objects.get(
            Q(cas_number__exact=cas_no) | Q(additional_cas__contains=cas_no)
        )
        data['object_exists'] = obj.pk
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
