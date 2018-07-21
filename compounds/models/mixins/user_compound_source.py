from bs4 import BeautifulSoup
from django.db import models
from django.core.exceptions import ValidationError
from django.contrib.postgres.fields import JSONField
from django.core.validators import RegexValidator


class UserCompoundSourceMixin(models.Model):
    """
    A mixin to provide fields which are common across various compound models and methods which set field data
    """

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
        validators=[RegexValidator(r'(?<![A-Z])[A-Z]{3}(?![A-Z])', "Format must be e.g. USD")],
    )
    amount = models.FloatField(
        max_length=20,
        help_text='Product amount in kg',
    )
    specification = models.CharField(
        max_length=100,
    )
    supplier = models.CharField(
        max_length=50,
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

    class Meta:
        abstract = True
