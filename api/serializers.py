from rest_framework import serializers
from .models import Compound, OdorType, Occurrence


class OdorTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = OdorType
        fields = ('term', )


class OccurrenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Occurrence
        fields = ('source', )


class CompoundSerializer(serializers.ModelSerializer):
    odor_categories = OdorTypeSerializer(many=True, read_only=True)
    occurrence = OccurrenceSerializer(read_only=True)

    class Meta:
        model = Compound
        fields = ('cas_number', 'smiles', 'iupac_name', 'trade_name', 'occurrence', 'odor_categories', )
