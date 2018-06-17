from rest_framework import serializers
from .models import Compound, OdorType


class OdorTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = OdorType
        fields = ('term', )


class CompoundSerializer(serializers.ModelSerializer):
    odor_categories = OdorTypeSerializer(many=True, read_only=True)

    class Meta:
        model = Compound
        fields = ('cas_number', 'smiles', 'iupac_name', 'trade_name', 'odor_categories', )
