from django.contrib.auth.mixins import LoginRequiredMixin
from django.http import HttpResponseForbidden
from django.views.generic import UpdateView

from compounds.models import Odorant
from compounds.forms import OdorantCompoundForm


class OdorantUpdateView(LoginRequiredMixin, UpdateView):
    model = Odorant
    form_class = OdorantCompoundForm
    template_name = 'odorants/odorant_update.html'
    context_object_name = 'compound'

    def get_initial(self):
        cpd = self.get_object()
        return {
            'iupac_name': cpd.iupac_name,
            'cid_number': cpd.cid_number,
            'chemical_name': cpd.chemical_name,
        }

    def post(self, request, *args, **kwargs):
        self.object = self.get_object()
        if not request.user.is_superuser and request.user.profile != self.object.edited_by:
            return HttpResponseForbidden()
        return super(OdorantUpdateView, self).post(request, *args, **kwargs)
