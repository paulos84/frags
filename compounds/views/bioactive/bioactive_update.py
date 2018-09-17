from django.contrib.auth.decorators import login_required
from django.http import HttpResponseForbidden
from django.utils.decorators import method_decorator
from django.views.generic import UpdateView

from compounds.forms import OdorantCompoundForm
from compounds.models import Odorant


class OdorantUpdateView(UpdateView):
    model = Odorant
    form_class = OdorantCompoundForm
    template_name = 'odorants/odorant_update.html'
    context_object_name = 'compound'

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_initial(self):
        cpd = self.get_object()
        return {
            'iupac_name': cpd.iupac_name,
            'cid_number': cpd.cid_number,
        }

    def post(self, request, *args, **kwargs):
        self.object = self.get_object()
        if not request.user.is_superuser or not request.user.profile == self.object.edited_by:
            return HttpResponseForbidden()
        return super(OdorantUpdateView, self).post(request, *args, **kwargs)
