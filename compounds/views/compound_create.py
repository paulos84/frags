from django.views.generic.edit import CreateView
from django.http import JsonResponse
from django.core.exceptions import ObjectDoesNotExist

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

def process_cas(request):
    username = request.GET.get('username', None)

    data = {
        'is_taken': Compound.objects.filter(username__iexact=username).exists()
    }
    try:
        cpd = Compound.objects.get(additional_cas__contains='177772-08-6')
        # redirect to corresponding detail view
    except ObjectDoesNotExist:
        pass

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