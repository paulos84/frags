from django.contrib.auth.mixins import LoginRequiredMixin
from django.shortcuts import redirect, reverse
from django.views.generic import ListView

from compounds.models import UserBioactive, UserOdorant


class UserAnalyticsView(LoginRequiredMixin, ListView):
    model = UserBioactive
    template_name = 'user/activity_list.html'
    context_object_name = 'user_bioactive_list'
    user_profile = None

    def dispatch(self, request, *args, **kwargs):
        if not request.user.is_superuser:
            return redirect(reverse('index'))
        return super(UserAnalyticsView, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super(UserAnalyticsView, self).get_context_data(**kwargs)
        return context
