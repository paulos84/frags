from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator

from django.views.generic import TemplateView, UpdateView, ListView
from django.views.generic.detail import SingleObjectMixin
from django.views.generic.edit import FormMixin
from compounds.forms import OdorantSearchForm, UserSourceCreateForm
from compounds.forms import UserLiteratureRefsForm
from compounds.models import Odorant, UserCompound, UserSource
from compounds.utils.find_literature import FindLiterature


# class UserSourcesView(SingleObjectMixin, TemplateView):
class UserSourceListView(LoginRequiredMixin, FormMixin, ListView):
    template_name = 'user/user_sources.html'
    user_compound = None
    context_object_name = 'sources_list'
    form_class = UserSourceCreateForm

    def dispatch(self, request, *args, **kwargs):
        try:
            self.user_compound = UserCompound.objects.get(
                user=request.user.profile,
                compound_id=kwargs['pk']
            )
        except UserCompound.DoesNotExist:
            pass
        return super(UserSourceListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        qs = UserSource.objects.filter(compound=self.user_compound)
        return qs

    def get_context_data(self, **kwargs):
        context = super(UserSourceListView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        context['form'] = self.form_class()
        return context

    def get_initial(self):
        initial = super(UserSourceListView, self).get_initial()
        initial['compound'] = self.user_compound
        return initial

    def post(self, request, *args, **kwargs):
        if request.POST:
            form = self.get_form()
            if form.is_valid():
                print(form.cleaned_data)
            else:
                print(form.errors)
