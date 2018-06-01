from django.views.generic.edit import UpdateView

from compounds.models import Compound
from compounds.forms import CompoundUpdateForm


class CompoundUpdateView(UpdateView):
    model = Compound
    form_class = CompoundUpdateForm
    template_name = 'compounds/update_compound.html'

    def get_context_data(self, **kwargs):
        context = super(CompoundUpdateView, self).get_context_data(**kwargs)
        context['iupac_name'] = self.get_object().iupac_name
        context['structure_url'] = self.get_object().structure_url
        return context