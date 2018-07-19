from django.shortcuts import reverse, HttpResponseRedirect
from django.views.generic import DeleteView

from compounds.models import UserCompound
from compounds.views.odorant.odorant_list import BaseOdorantListView


# TODO if notes go beyond area use ... truncation


class ActivityListView(BaseOdorantListView):
    pass


class UserCompoundNotesDeleteView(DeleteView):
    model = UserCompound

    def delete(self, request, *args, **kwargs):
        self.object = self.get_object()
        if self.request.user.is_superuser or self.request.user.profile == self.object.user:
            self.object.delete()
            # messages.success(self.request, self.success_message)
        return HttpResponseRedirect(self.get_success_url())

    def get_success_url(self):
        return reverse(
            'odorant-detail', kwargs={'pk': self.object.compound.pk})

