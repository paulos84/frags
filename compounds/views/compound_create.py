from django.urls import reverse_lazy
from django.views.generic.edit import CreateView, DeleteView, UpdateView

from compounds.models.compound import Compound

# Bottom of page: using AJAX with formviews:
## https://docs.djangoproject.com/en/2.0/topics/class-based-views/generic-editing/

class CompoundUpdate(UpdateView):
    model = Compound
    fields = ['name']

