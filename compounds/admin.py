from django.contrib import admin

from compounds.models.compound import Compound
from compounds.models.odor import Odor


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    pass


@admin.register(Odor)
class OdorAdmin(admin.ModelAdmin):
    pass
