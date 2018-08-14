from django.template import Library

register = Library()


@register.filter
def user_notes(bioactive, user_profile):
    return bioactive.user_activities(user_profile)['notes']


@register.filter
def user_sources(bioactive, user_profile):
    return bioactive.user_activities(user_profile)['sources']


@register.filter
def user_lit_refs(bioactive, user_profile):
    return bioactive.user_activities(user_profile)['lit_refs']
