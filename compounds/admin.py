from django.contrib import admin

from compounds.models import Compound, Occurrence, OdorType


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    pass


@admin.register(Occurrence)
class OccurrenceAdmin(admin.ModelAdmin):
    pass


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass
