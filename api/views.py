from rest_framework.views import APIView
from rest_framework.response import Response
from .serializers import BioactiveSerializer, CompoundSerializer
from .models import Odorant, Bioactive
from rest_framework_swagger.views import get_swagger_view

schema_view = get_swagger_view(title='Fragrance odorants and materials API', url=None)


class AllCompoundsData(APIView):
    def get(self, request, format=None):
        queryset = Odorant.objects.all()
        serializer = CompoundSerializer(queryset, many=True)
        return Response(serializer.data)


class AllBioactivesData(APIView):
    def get(self, request, format=None):
        queryset = Bioactive.objects.all()
        serializer = BioactiveSerializer(queryset, many=True)
        return Response(serializer.data)
