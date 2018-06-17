from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver

from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType

from compounds.models import Compound


class Profile(models.Model):
    """ How to use   https://simpleisbetterthancomplex.com/tutorial/2016/07/22/how-to-extend-django-user-model.html  """

    user = models.OneToOneField(User, on_delete=models.CASCADE)
    bio = models.TextField(max_length=500, blank=True)
    location = models.CharField(max_length=30, blank=True)
    birth_date = models.DateField(null=True, blank=True)


@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        Profile.objects.create(user=instance)


@receiver(post_save, sender=User)
def save_user_profile(sender, instance, **kwargs):
    instance.profile.save()


class CompoundNotes(models.Model):
    notes = models.TextField(
        max_length=500,
    )
    compound = models.ForeignKey(
        Compound,
        on_delete=models.CASCADE,
    )
    user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
    )
    # activities = GenericRelation(Activity)


# class Activity(models.Model):
#     ACTIVITY_TYPES = (
#         (1, 'Compound notes'),
#     )
#     user = models.ForeignKey(
#         Profile,
#         on_delete=models.CASCADE
#     )
#     activity_type = models.IntegerField(choices=ACTIVITY_TYPES)
#
#     content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
#     object_id = models.PositiveIntegerField()
#     content_object = GenericForeignKey()
