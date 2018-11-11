from django.contrib import admin

from compounds.forms.admin import ActivityAdminForm, BioactiveAdminForm, SubstructureAdminForm
from compounds.models import (Activity, Bioactive, BioactiveCore, Odorant, UserOdorant, OdorType, Profile, Substructure,
                              Enzyme, CompoundSource, UserBioactive)


@admin.register(Activity)
class ActivityAdmin(admin.ModelAdmin):
    form = ActivityAdminForm

    def get_queryset(self, request):
        qs = super(ActivityAdmin, self).get_queryset(request)
        return qs.order_by('action__name', 'name')


@admin.register(Bioactive)
class BioactiveAdmin(admin.ModelAdmin):
    form = BioactiveAdminForm

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('chemical_properties', 'iupac_name', 'smiles')
        return self.readonly_fields


@admin.register(BioactiveCore)
class BioactiveCoreAdmin(admin.ModelAdmin):

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('cid_number', 'iupac_name', )
        return self.readonly_fields


@admin.register(Odorant)
class OdorantAdmin(admin.ModelAdmin):

    def get_readonly_fields(self, request, obj=None):
        if obj:  # editing an existing object
            return self.readonly_fields + ('cas_number', 'chemical_properties', 'iupac_name', 'smiles')
        return self.readonly_fields


@admin.register(OdorType)
class OdorAdmin(admin.ModelAdmin):
    pass


@admin.register(Enzyme)
class EnzymeAdmin(admin.ModelAdmin):
    pass


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    pass


@admin.register(Substructure)
class SubstructureAdmin(admin.ModelAdmin):
    form = SubstructureAdminForm
    readonly_fields = ['iupac_name']


@admin.register(UserBioactive)
class UserBioactiveAdmin(admin.ModelAdmin):
    pass


@admin.register(UserOdorant)
class UserOdorantAdmin(admin.ModelAdmin):
    pass


@admin.register(CompoundSource)
class UserOdorantSourceAdmin(admin.ModelAdmin):
    pass
