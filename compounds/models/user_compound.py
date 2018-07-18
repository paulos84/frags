from django.db import models
from django.contrib.postgres.fields import ArrayField

from compounds.models import Odorant, Profile
from .managers.activity import ActivityManager


class UserCompound(models.Model):
    notes = models.TextField(
        max_length=500,
        null=True,
        blank=True,
    )
    compound = models.ForeignKey(
        Odorant,
        on_delete=models.CASCADE,
        related_name='notes_set',
        blank=True,
    )
    user = models.ForeignKey(
        Profile,
        on_delete=models.CASCADE,
        related_name='notes_set',
        blank=True,
    )
    literature_refs = ArrayField(
        models.CharField(max_length=100),
        null=True,
        blank=True,
    )

    objects = ActivityManager()

    def __str__(self):
        return 'notes: ' + str(self.compound) + '_' + str(self.user)

    @classmethod
    def lit_refs_actions(cls, request, refs, compound):
        message = None
        if 'save_refs' in request.POST:
            instance, _ = cls.objects.get_or_create(
                user=request.user.profile,
                compound=compound
            )
            if instance.literature_refs:
                instance.literature_refs.extend(refs)
            else:
                instance.literature_refs = refs
            message = 'Article{} saved'.format('s' if len(refs) > 1 else '')
            instance.save()
        if 'remove_refs' in request.POST:
            instance = UserCompound.objects.get(
                user=request.user.profile,
                compound=compound
            )
            for ref in refs:
                instance.literature_refs.remove(ref)
            message = 'Article{} removed'.format('s' if len(refs) > 1 else '')
            instance.save()
        return message

