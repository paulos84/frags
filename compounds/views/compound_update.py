from django.shortcuts import redirect
from django.views.generic import UpdateView
from django.utils import timezone

from compounds.models import Compound
from compounds.forms import CompoundUpdateForm


class CompoundUpdateView(UpdateView):
    model = Compound
    form_class = CompoundUpdateForm
    template_name = 'compounds/compound_update.html'
    # success_url = "/accounts/profile/"
    context_object_name = 'post'

    # def form_valid(self, form):
    #     if form.is_valid:
    #         user = form.save()
    #         user = authenticate(username=user.username, password=form.cleaned_data['password1'])
    #         login(self.request, user)
    #         return super(RegistrationView, self).form_valid(form)  # I still have no idea what this is
