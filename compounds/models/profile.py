from django.contrib.auth.models import User
from django.db import models
from django.db.models.signals import post_save
from django.dispatch import receiver


class Profile(models.Model):
    """
    A model representing a User and their activity
    """
    user = models.OneToOneField(
        User,
        on_delete=models.CASCADE,
    )
    email_confirmed = models.BooleanField(
        default=False,
    )
    bioactive_set = models.ManyToManyField(
        'compounds.Bioactive',
        through='compounds.UserBioactive',
        verbose_name='Bioactive compound set',
        blank=True,
    )
    odorant_set = models.ManyToManyField(
        'compounds.Odorant',
        through='compounds.UserOdorant',
        verbose_name='Bioactive compound set',
        blank=True,
    )

    def __str__(self):
        return self.user.username + '_profile'


@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
    """
    Creates an instance when a User instance is created
    """
    if created:
        Profile.objects.create(user=instance)


@receiver(post_save, sender=User)
def save_user_profile(sender, instance, **kwargs):
    """
    Updates an instance when a User instance is updated
    """
    instance.profile.save()
