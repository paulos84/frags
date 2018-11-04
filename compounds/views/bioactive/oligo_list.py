from django.contrib import messages
from django.shortcuts import redirect, reverse, get_object_or_404
from django.utils.text import slugify
from django.views.generic import ListView

from compounds.models import Activity, Bioactive, BioactiveCore
from compounds.views.mixins import BioactiveSearchFilterMixin


class OligosaccharideListView(BioactiveSearchFilterMixin, ListView):
    template_name = 'bioactives/bioactive_list.html'
    context_object_name = 'bioactive_list'
    paginate_by = 32

    def get_queryset(self):
        biocore = get_object_or_404(BioactiveCore, name='Oligosaccharides')
        return biocore.bioactives.all()

    def get_context_data(self, **kwargs):
        context = super(OligosaccharideListView, self).get_context_data(**kwargs)
        context.update({
            'page_header': 'Oligosaccharides',
            'category': 2,
        })
        return context



# except Bioactive.DoesNotExist:
# messages.info(request, 'No compound matching InChiKey')
# return redirect(reverse('bioactive-list', kwargs={'category': 'medicinal'}))
