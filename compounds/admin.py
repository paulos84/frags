from django.contrib import admin

from compounds.models.compound import Compound
from compounds.models.odor_type import OdorType


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    pass


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass
