from django.db import models
from django.contrib.postgres.fields import ArrayField


class UserCompoundMixin(models.Model):
    notes = models.TextField(
        max_length=500,
        null=True,
        blank=True,
    )
    literature_refs = ArrayField(
        models.CharField(max_length=100),
        null=True,
        blank=True,
    )

    class Meta:
        abstract = True

    @classmethod
    def lit_refs_actions(cls, request, refs, compound):
        if 'save_refs' in request.POST:
            instance, _ = cls.objects.get_or_create(
                user=request.user.profile,
                compound=compound
            )
            if instance.literature_refs:
                instance.literature_refs.extend(refs)
            else:
                instance.literature_refs = refs
            instance.save()
        elif 'remove_refs' in request.POST:
            instance = cls.objects.get(
                user=request.user.profile,
                compound=compound
            )
            for ref in refs:
                instance.literature_refs.remove(ref)
            instance.save()
