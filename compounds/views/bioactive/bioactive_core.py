from django.db.models import Q
from django.views.generic import ListView
from django.views.generic.detail import SingleObjectMixin

from compounds.models import Odorant, UserOdorant, OdorType, Substructure
from compounds.forms import OdorantSearchForm
from compounds.views.mixins.search_filter import OdorantSearchFilterMixin


class BioactiveCoreMatchList(OdorantSearchFilterMixin, SingleObjectMixin, ListView):
    paginate_by = 14
    template_name = "odorants/substructure_detail.html"

    def get(self, request, *args, **kwargs):
        self.object = self.get_object(queryset=Substructure.objects.all())
        return super(BioactiveCoreMatchList, self).get(request, *args, **kwargs)

    def get_queryset(self):
        return self.object.odorant_set()