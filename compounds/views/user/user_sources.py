from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import get_object_or_404, redirect
from django.urls import reverse
from django.utils.decorators import method_decorator

from django.views.generic import TemplateView, UpdateView, ListView
from django.views.generic.detail import SingleObjectMixin
from django.views.generic.edit import FormMixin
from compounds.forms import OdorantSearchForm, UserOdorantSourceCreateForm
from compounds.forms import UserLiteratureRefsForm
from compounds.models import Odorant, UserOdorant, UserOdorantSource
from compounds.utils.find_literature import FindLiterature


class UserOdorantSourceListView(LoginRequiredMixin, FormMixin, ListView):
    template_name = 'user/user_sources.html'
    user_compound = None
    odorant = None
    context_object_name = 'sources_list'
    form_class = UserOdorantSourceCreateForm

    def dispatch(self, request, *args, **kwargs):
        self.odorant = get_object_or_404(Odorant, pk=kwargs['pk'])
        try:
            self.user_compound = UserOdorant.objects.get(
                user=request.user.profile,
                compound_id=kwargs['pk']
            )
        except UserOdorant.DoesNotExist:
            pass
        return super(UserOdorantSourceListView, self).dispatch(request, *args, **kwargs)

    def get_queryset(self):
        qs = UserOdorantSource.objects.filter(
            compound=self.user_compound,
            user=self.request.user.profile)
        return qs

    def get_context_data(self, **kwargs):
        context = super(UserOdorantSourceListView, self).get_context_data(**kwargs)
        context['compound_search'] = OdorantSearchForm()
        context['form'] = self.form_class()
        return context

    def post(self, request, *args, **kwargs):
        if request.POST:
            form = self.get_form()
            if form.is_valid():
                print(form.cleaned_data)
                return self.form_valid(form)
            else:
                print(form.errors)
                return self.form_invalid(**kwargs)

    def form_valid(self, form):
        form.instance.odorant = self.odorant
        form.instance.user = self.request.user.profile
        form.instance.save()
        return super(UserOdorantSourceListView, self).form_valid(form)

    def get_success_url(self):
        return reverse('user-odorant-sources', kwargs={'pk': self.kwargs['pk']})
