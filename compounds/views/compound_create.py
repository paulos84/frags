from django.views.generic.edit import CreateView
from django.http import JsonResponse, HttpResponseRedirect
from django.shortcuts import redirect
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Q
import cirpy
import pubchempy as pcp

from compounds.models.compound import Compound
from compounds.forms.compound_forms import CompoundCreateForm


class CompoundCreateView(CreateView):
    model = Compound
    form_class = CompoundCreateForm
    template_name = 'compounds/create_compound.html'

    # def get_context_data(self, **kwargs):
    #     context = super(CompoundCreateForm, self).get_context_data(**kwargs)
    #     cid_no = self.get_object().cid_number
    #     context['structure_url'] = 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no)
    #     try:
    #         context['synonyms'] = ', '.join(pcp.get_compounds(cid_no)[0].synonyms)
    #         # PARSE OUT EC... from synonyms and FEMA...
    #     except KeyError:
    #         context['synonyms'] = 'n/a'
    #     return context

# Todo: try and keep data from the lookup below on the object in form so model doesn't have to do api calls again

def process_cas(request):
    cas_no = request.GET.get('cas_number')
    # print('foo')
    # data = {}
    # try:
    #     obj = Compound.objects.get(
    #         Q(cas_number__exact=cas_no) | Q(additional_cas__contains=cas_no)
    #     )
    #     return redirect(obj)
    # except ObjectDoesNotExist:
    #     pass
    # try:
    #     smiles = cirpy.query(cas_no, 'smiles')[0].value
    #     data['iupac_name'] = cirpy.Molecule(smiles).iupac_name
    #     cid_no = pcp.get_compounds(smiles, 'smiles')[0].cid
    #     data['structure_url'] = 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no)
    # except IndexError:
    #     data['error_message'] = 'A user with this username already exists.'
    data = {
        'is_taken': True
    }
    return JsonResponse(data)



    # cid_no = self.get_object().cid_number
    # context['structure_url'] = 'https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={}&amp;t=l'.format(cid_no)


    #
    # def validate_username(request):
    #     username = request.GET.get('username', None)
    #     data = {
    #         'is_taken': User.objects.filter(username__iexact=username).exists()
    #     }
    #     if data['is_taken']:
    #         data['error_message'] = 'A user with this username already exists.'
    #     return JsonResponse(data)




#ModelForm.save for cleaning form data, Model.save for cleaning object attributes.

# These views inherit SingleObjectTemplateResponseMixin which uses template_name_suffix to construct the template_name based on the model.
#
# In this example:
#
#     CreateView and UpdateView use myapp/author_form.html
#     DeleteView uses myapp/author_confirm_delete.html
#
# If you wish to have separate templates for CreateView and UpdateView, you can set either template_name or template_name_suffix on your v