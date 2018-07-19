from django.contrib.auth.mixins import LoginRequiredMixin
from django.db.models import Q

from compounds.models import Bioactive, UserCompound
from compounds.views.base_compound_list import BaseCompoundListView


class BaseBioactiveListView(BaseCompoundListView):
    queryset = Bioactive.objects.all()
    template_name = 'odorants/odorant_list.html'
    paginate_by = 32

    def get_context_data(self, **kwargs):
        context = super(BaseBioactiveListView, self).get_context_data(**kwargs)
        # context['odor_types'] = OdorType.objects.values('term')
        # some other category...e.g. func food, medicinal,
        return context


class BioactiveListView(BaseBioactiveListView):
    usage_type = None

    def dispatch(self, request, *args, **kwargs):
        self.usage_type = kwargs['usage']
        return super(BioactiveListView, self).dispatch(request, *args, **kwargs)

    # SET A USAGE_TYPE CHOICES INTFIELD ON MODEL AND USE THIS PROPERTY TO FILTER (STILL USE DISPATCH/KWARG?)

    def get_context_data(self, **kwargs):
        context = super(BioactiveListView, self).get_context_data(**kwargs)
        # context['odor_types'] = OdorType.objects.values('term')
        # some other category...e.g. func food, medicinal,
        return context

    def get_queryset(self):
        queryset = super(BioactiveListView, self).get_queryset()
        # filter based upon usage_type set with kwargs in dispatch
        return queryset


class UserBioactiveListView(LoginRequiredMixin, BaseBioactiveListView):
    template_name = 'odorants/user_odorant_list.html'
    context_object_name = 'compound_list'
    usage_type = None

    def get_queryset(self):
        queryset = super(UserBioactiveListView, self).get_queryset()
        # filter based upon usage_type set with kwargs in dispatch
        return queryset

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['page_header'] = 'My compound notes'
        return context
