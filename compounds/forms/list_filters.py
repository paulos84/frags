from django.forms import forms
import django_filters

from compounds.models import Compound, OdorType


class CompoundFilter(django_filters.FilterSet):
    iupac_name = django_filters.CharFilter(
        lookup_expr='icontains', label='IUPAC name contains', )
    trade_name = django_filters.CharFilter(
        lookup_expr='iexact', )
    odor_description = django_filters.CharFilter(
        lookup_expr='icontains', label='Scent keywords', )

    class Meta:
        model = Compound
        fields = [
            'iupac_name', 'trade_name', 'odor_description',
        ]