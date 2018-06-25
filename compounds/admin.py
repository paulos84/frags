from django.contrib import admin

from compounds.models import Compound, UserNotes, OdorType, Profile, Substructure


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    pass


@admin.register(UserNotes)
class CompoundNotesAdmin(admin.ModelAdmin):
    pass


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    pass


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass


@admin.register(Substructure)
class SubstructureAdmin(admin.ModelAdmin):
    pass
