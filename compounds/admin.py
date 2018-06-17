from django.contrib import admin

from compounds.models import Compound, OdorType, Substructure


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    pass


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass


@admin.register(Substructure)
class SubstructureAdmin(admin.ModelAdmin):
    pass
