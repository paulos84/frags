from django.db import models
from django.core.validators import RegexValidator
from django.core.exceptions import ValidationError


class CompoundSource(models.Model):

    user_odorant = models.ForeignKey(
        'compounds.UserOdorant',
        related_name='userodorant_sources',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    user_bioactive = models.ForeignKey(
        'compounds.UserBioactive',
        related_name='userbioactive_sources',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    price = models.DecimalField(
        decimal_places=2,
        max_digits=10,
    )
    currency_choices = (
        ('', ''),
        ('USD', 'US Dollars'),
        ('HKD', 'Hong Kong Dollar'),
        ('CNY', 'Chinese Yuan'),
        ('JPY', 'Japanese Yen'),
        ('GBP', 'British Pound'),
        ('EUR', 'Euro'),
        ('IND', 'Indian Rupee'),
    )
    currency = models.CharField(
        max_length=3,
        choices=currency_choices,
        validators=[RegexValidator(r'(?<![A-Z])[A-Z]{3}(?![A-Z])', "Format must be e.g. USD")],
        blank=True,
        default='',
    )
    amount = models.FloatField(
        max_length=20,
        help_text='Product amount in g',
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

    def __str__(self):
        return '{currency} {price} / {amount}{unit}'.format(
            currency=self.currency, price=self.price, amount=self.amount, unit='kg')

    def save(self, *args, **kwargs):
        if not any([self.user_odorant, self.user_bioactive]):
            raise ValidationError('Something went wrong B')
        super(CompoundSource, self).save(*args, **kwargs)
