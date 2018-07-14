from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import viewsets
from .serializers import CompoundSerializer
from datetime import datetime, date, timedelta
from .models import Odorant
from rest_framework_swagger.views import get_swagger_view

schema_view = get_swagger_view(title='Fragrance odorants and materials API', url=None)


class AllCompoundsData(APIView):

    def get(self, request, format=None):
        """ filter results according to the site code """
        queryset = Odorant.objects.all()
        serializer = CompoundSerializer(queryset, many=True)
        return Response(serializer.data)


# class RecentSiteData(APIView):
#
#     def get(self, request, code, days, format=None):
#         """ filter results according to the site code and number of recent days """
#         queryset = Compound.objects.filter(site__code=code)
#         serializer = CompoundSerializer(queryset, many=True)
#         return Response(serializer.data)
