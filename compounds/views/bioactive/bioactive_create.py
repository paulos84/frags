from django.views.generic.edit import CreateView
from django.http import JsonResponse
from django.core.exceptions import ObjectDoesNotExist
from django.shortcuts import redirect
import cirpy
import pubchempy as pcp

from compounds.models import Activity, Bioactive
from compounds.forms import BioactiveCreateForm, OdorantSearchForm


class BioactiveCreateView(CreateView):
    model = Bioactive
    form_class = BioactiveCreateForm
    template_name = 'bioactives/bioactive_create.html'

    def get_context_data(self, **kwargs):
        context = super(BioactiveCreateView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        return context

    def form_valid(self, form):
        return super().form_valid(form)


def process_bioactive_identifier(request):
    cas_no = request.GET.get('cas_number')
    inchikey = request.GET.get('inchikey')
    if cas_no:
        obj = Bioactive.objects.filter(chemical_properties__synonyms__icontains=cas_no).first()
    elif inchikey:
        obj = Bioactive.objects.filter(inchikey__exact=inchikey).first()
    if obj:
        data = {
            'object_exists': obj.get_absolute_url(),
            'object_exists_name': obj.chemical_name,
        }
        return JsonResponse(data)
    try:
        if cas_no:
            smiles = cirpy.query(cas_no, 'smiles')[0].value
            pcp_query = pcp.get_compounds(smiles, 'smiles')[0]
            if not pcp_query.cid:
                raise IndexError
        else:
            pcp_query = pcp.get_compounds(inchikey, 'inchikey')[0]
            if not pcp_query.cid:
                raise IndexError
    except (IndexError, pcp.BadRequestError):
        return JsonResponse({
            'error': 'No compound found for this CAS number'
        })
    data = {
        'chemical_name': Bioactive.scrape_compound_name(pcp_query.cid),
        'iupac_name': pcp_query.iupac_name,
        'inchikey': pcp_query.inchikey,
        'structure_url': 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(pcp_query.cid),
        'hidden_cid': pcp_query.cid,
        'smiles': pcp_query.isomeric_smiles or pcp_query.canonical_smiles or '',
    }
    return JsonResponse(data)


def process_activity(request):
    print(request.GET)
    category_choice = request.GET.get('category_choice')
    classification_choice = request.GET.get('classification_1')
    action_choice = request.GET.get('action')
    if category_choice:
        applicable_classifications = Activity.applicable_classifications(category_choice)
        categories = [{'name': a[1]} for a in applicable_classifications]
        categories.insert(0, {'name': '-------'})
        return JsonResponse({'categories': categories,'has_children': True}, safe=False)
    if classification_choice:
        choice = Activity.map_to_classification(classification_choice)
        activities = list(Activity.objects.actions().filter(
            classification=choice).values('name'))
        activities.insert(0, {'name': '-------'})
        return JsonResponse({'actions': activities, 'has_children': True}, safe=False)
    if action_choice:
        selected_action = Activity.objects.get(
            category=1,
            classification=Activity.classifications[int(request.GET.get('parent_classification')) - 1][0]
        )
        mechanisms = list(selected_action.mechanisms.values('name'))
        print(mechanisms)
        return JsonResponse(mechanisms, safe=False)




    # TODO: check if has mechanisms...is has..unhide mechanisms field






# def process_cas(request):
#     data = {}
#     cas_no = request.GET.get('cas_number')
#     try:
#         obj = Odorant.objects.get(cas_number__exact=cas_no)
#         data['object_exists'] = obj.get_absolute_url()
#         data['object_exists_name'] = obj.iupac_name
#         return JsonResponse(data)
#     except ObjectDoesNotExist:
#         pass
#     try:
#         smiles = cirpy.query(cas_no, 'smiles')[0].value
#         cid_no = pcp.get_compounds(smiles, 'smiles')[0].cid
#         if smiles and cid_no:
#             data['iupac_name'] = cirpy.Molecule(smiles).iupac_name
#             data['structure_url'] = 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no)
#             data.update({'smiles': smiles, 'hidden_cid': cid_no})
#     except IndexError:
#         data['error'] = 'No compound found for this CAS number'
#     return JsonResponse(data)
