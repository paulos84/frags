from django.urls import reverse_lazy
from django.views.generic.edit import UpdateView

from compounds.models import Compound
from compounds.forms.compound_forms import CompoundUpdateForm


# Bottom of page: using AJAX with formviews:
## https://docs.djangoproject.com/en/2.0/topics/class-based-views/generic-editing/


class CompoundUpdateView(UpdateView):
    model = Compound
    form_class = CompoundUpdateForm
    template_name = 'compounds/update_compound.html'





#ModelForm.save for cleaning form data, Model.save for cleaning object attributes.

# These views inherit SingleObjectTemplateResponseMixin which uses template_name_suffix to construct the template_name based on the model.
#
# In this example:
#
#     CreateView and UpdateView use myapp/author_form.html
#     DeleteView uses myapp/author_confirm_delete.html
#
# If you wish to have separate templates for CreateView and UpdateView, you can set either template_name or template_name_suffix on your view class.