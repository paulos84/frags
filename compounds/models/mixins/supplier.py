from django.db import models
from django.core.exceptions import ValidationError


class SupplierMixin(models.Model):

    choices = (
        ('BASF', 'BASF AG, Germany'),
        ('Danisco', 'Danisco A/S, Denmark'),
        ('Firmenich', 'Firmenich SA, Switzerland'),
        ('Giv.', 'Givaudan SA, Switzerland'),
        ('IFF', 'International Flavors & Fragrances, USA'),
        ('Quest', 'Quest International, UK'),
        ('Symrise', 'Symrise GmbH & Co KG, Germany'),
        ('Takasago', 'Takasago Perfumery Co., Japan'),
        ('Vioryl', 'Vioryl SA, Greece'),
    )
    supplier = models.CharField(
        max_length=25, default='',
        choices=choices,
        blank=True,
    )
    trade_name = models.CharField(
        max_length=20, default='',
        verbose_name='Trade name',
        blank=True,
    )

    class Meta:
        abstract = True

    def clean(self):
        if self.supplier is not None and self.trade_name is None:
            raise ValidationError('Trade name required if a supplier is entered')
        super(SupplierMixin, self).clean()
