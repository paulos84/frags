from django.db import models


# CSV import of .. price..amount..specification..weblink?

class UserSource(models.Model):

    price = models.DecimalField(
        max_length=10,
    )
    currency = 
    amount = models.FloatField(
        max_length=20,
        help_text='Product amount in kg'
    )
    specification = models.CharField(
        max_length=100,
        default='',
        blank=True,
    )
    supplier = models.CharField(
        max_length=50,
        default='',
        blank=True,
    )
    product_number = models.CharField(
        max_length=20,
        default='',
        blank=True,
    )
    url = models.URLField(
        max_length=500,
        default='',
        blank=True,
    )
    compound = models.ManyToManyField(
        'compounds.UserCompound',
        related_name='sources',
        verbose_name='User compounds sources',
    )

