from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import reverse, HttpResponseRedirect
from django.views.generic import DeleteView

from compounds.models import UserBioactive, UserOdorant
from compounds.views.odorant.odorant_list import BaseOdorantListView


class UserActivityListView(LoginRequiredMixin, BaseOdorantListView):
    model = UserBioactive
    template_name = 'user/activity_list.html'
    context_object_name = 'bioactive_list'
    user_profile = None

    def dispatch(self, request, *args, **kwargs):
        self.user_profile = self.request.user.profile
        self.queryset = self.user_profile.bioactive_set.all()
        return super(UserActivityListView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(UserActivityListView, self).get_context_data(**kwargs)
        context['odorant_list'] = self.user_profile.odorant_set.all()
        return context


class UserCompoundNotesDeleteView(DeleteView):
    model = None

    def dispatch(self, request, *args, **kwargs):
        model_kwarg = kwargs.pop('model')
        self.model = UserOdorant if model_kwarg == 'odorant' else UserBioactive
        return super(UserCompoundNotesDeleteView, self).dispatch(request, *args, **kwargs)

    def delete(self, request, *args, **kwargs):
        self.object = self.get_object()
        if self.request.user.is_superuser or self.request.user.profile == self.object.user:
            self.object.notes = ''
            self.object.save()
        return HttpResponseRedirect(self.get_success_url())

    def get_success_url(self):
        name = 'odorant-detail' if self.model == UserOdorant else 'bioactive-detail'
        return reverse(
            name, kwargs={'pk': self.object.compound.pk}
        )
