import django_filters

from compounds.models import Compound


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

    def __init__(self, *args, **kwargs):
        super(CompoundFilter, self).__init__(*args, **kwargs)
        # at sturtup user doen't push Submit button, and QueryDict (in data) is empty
        if self.data == {}:
            self.queryset = self.queryset.none()