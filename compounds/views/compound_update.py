from django.views.generic import UpdateView
from django.utils import timezone
from django.contrib.auth.decorators import login_required
from django.utils.decorators import method_decorator
from django.http import HttpResponseForbidden

from compounds.models import Compound
from compounds.forms import CompoundUpdateForm, EditCompoundForm


class CompoundUpdateView(UpdateView):
    model = Compound
    form_class = EditCompoundForm
    template_name = 'compounds/compound_update.html'
    # success_url = "/accounts/profile/"
    context_object_name = 'compound'

    @method_decorator(login_required)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    # def form_valid(self, form):
    #     if form.is_valid:
    #         user = form.save()
    #         user = authenticate(username=user.username, password=form.cleaned_data['password1'])
    #         login(self.request, user)
    #         return super(RegistrationView, self).form_valid(form)  # I still have no idea what this is

    def get_initial(self):
        cpd = self.get_object()
        return {
            'iupac_name': cpd.iupac_name,
            'cid_number': cpd.cid_number,
        }

    def post(self, request, *args, **kwargs):
        self.object = self.get_object()
        if not request.user.is_superuser or not request.user.profile == self.object.created_by:
            return HttpResponseForbidden()
        return super(CompoundUpdateView, self).post(request, *args, **kwargs)