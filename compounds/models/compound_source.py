from django.db import models
from django.core.validators import RegexValidator
from django.core.exceptions import ValidationError


class CompoundSource(models.Model):

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
    unit_choices = (
        ('g', 'grams'),
        ('kg', 'kilograms'),
    )
    unit = models.CharField(
        max_length=2,
        choices=unit_choices,
        blank=True,
        default='g',
    )
    amount = models.FloatField(
        max_length=10,
    )
    specification = models.CharField(
        max_length=50,
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
        max_length=100,
        default='',
        blank=True,
    )
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
    odorant = models.ForeignKey(
        'compounds.Odorant',
        related_name='odorant_sources',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    bioactive = models.ForeignKey(
        'compounds.Bioactive',
        related_name='bioactive_sources',
        on_delete=models.CASCADE,
        null=True,
        blank=True,
    )
    source = models.ForeignKey(
        'self',
        on_delete=models.CASCADE,
        blank=True,
        null=True,
    )

    def __str__(self):
        return '{} | {} | {}'.format(
            self.user_odorant or self.user_bioactive, self.supplier, self.amount)

    @property
    def summary_display(self):
        return '{currency} {price} / {amount}{unit}'.format(
            currency=self.currency, price=self.price, amount=self.amount, unit=self.unit)

    def save(self, *args, **kwargs):
        if not any([self.user_odorant, self.user_bioactive]):
            raise ValidationError('Something went wrong B')
        super(CompoundSource, self).save(*args, **kwargs)
