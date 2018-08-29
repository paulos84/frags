from rest_framework import serializers
from .models import Odorant, OdorType


class OdorTypeSerializer(serializers.ModelSerializer):
    class Meta:
        model = OdorType
        fields = ('term', )


class CompoundSerializer(serializers.ModelSerializer):
    odor_categories = OdorTypeSerializer(many=True, read_only=True)

    class Meta:
        model = Odorant
        fields = ('cas_number', 'smiles', 'iupac_name', 'chemical_name', 'chemical_properties', 'odor_categories', )


class BioactiveSerializer(serializers.ModelSerializer):
    bioactives = OdorTypeSerializer(many=True, read_only=True)

    class Meta:
        model = Odorant
        fields = ('inchikey', 'smiles', 'iupac_name', 'chemical_name', 'chemical_properties', )

