from django.contrib.auth.decorators import login_required
from django.shortcuts import HttpResponseRedirect
from django.urls import reverse
from django.utils.decorators import method_decorator
from django.views.generic import DetailView
from django.views.generic.edit import FormMixin

from compounds.models import Bioactive, UserCompound
from compounds.forms import CompoundNotesForm, OdorantUpdateForm
from compounds.views.base_compound_list import BaseCompoundListView


class BioactiveDetailView(DetailView):
    model = Bioactive
    template_name = 'bioactives/bioactive_detail.html'

    def get_context_data(self, **kwargs):
        context = super(BioactiveDetailView, self).get_context_data(**kwargs)
        return context

    def get_success_url(self):
        return reverse('bioactive-detail', kwargs={'pk': self.object.pk})
