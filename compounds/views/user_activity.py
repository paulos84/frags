from django.views.generic import DeleteView, UpdateView
from django.shortcuts import reverse, HttpResponseRedirect

from .compound_list import BaseCompoundListView
from compounds.models import Odorant, UserNotes

# TODO if notes go beyond area use ... truncation


class ActivityListView(BaseCompoundListView):
    pass


class UserCompoundNotesDeleteView(DeleteView):
    model = UserNotes

    def delete(self, request, *args, **kwargs):
        self.object = self.get_object()
        if self.request.user.is_superuser or self.request.user.profile == self.object.user:
            self.object.delete()
            # messages.success(self.request, self.success_message)
        return HttpResponseRedirect(self.get_success_url())

    def get_success_url(self):
        return reverse(
            'odorant-detail', kwargs={'pk': self.object.compound.pk})

