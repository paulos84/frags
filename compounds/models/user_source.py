from django.db import models
from django.core.validators import RegexValidator

# CSV import of .. price..amount..specification..weblink?

class UserSource(models.Model):

    currency_choices = (
        ('USD', 'US Dollars'),
        ('HKD', 'Hong Kong Dollar'),
        ('CNY', 'Chinese Yuan'),
        ('JPY', 'Japanese Yen'),
        ('GBP', 'British Pound'),
        ('EUR', 'Euro'),
        ('IND', 'Indian Rupee'),
    )

    price = models.DecimalField(
        decimal_places=2,
        max_digits=10,
    )
    currency = models.CharField(
        max_length=3,
        choices=currency_choices,
        help_text='Currency code',
        validators=[RegexValidator(r'\d+(?:-\d+)+', "Format must be e.g. USD")],
    )
    amount = models.FloatField(
        max_length=20,
        help_text='Product amount in kg',
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

