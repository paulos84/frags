from django.forms import forms
import django_filters

from compounds.models import Compound


class CompoundFilter(django_filters.FilterSet):
    iupac_name = django_filters.CharFilter(lookup_expr='icontains')
    trade_name = django_filters.CharFilter(lookup_expr='iexact')

    class Meta:
        model = Compound
        fields = ['cas_number', 'iupac_name', 'trade_name']