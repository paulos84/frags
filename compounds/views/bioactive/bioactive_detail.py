from django.urls import reverse
from django.views.generic import DetailView

from compounds.models import Bioactive


class BioactiveDetailView(DetailView):
    model = Bioactive
    template_name = 'bioactives/bioactive_detail.html'

    def get_context_data(self, **kwargs):
        context = super(BioactiveDetailView, self).get_context_data(**kwargs)
        return context

    def get_success_url(self):
        return reverse('bioactive-detail', kwargs={'pk': self.object.pk})
