from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError
from django.db.models import Q

from compounds.models import Compound


class CompoundUpdateForm(forms.ModelForm):

    class Meta:
        model = Compound
        fields = ['odor_description', 'odor_categories', 'supplier', 'trade_name']
        initial = {'odor_description': 'foo', 'trade_name': 'bar' }
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 5, 'cols': 52, }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '8', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
        }

    def __init__(self, *args, **kwargs):
        super(CompoundUpdateForm, self).__init__(*args, **kwargs)
        self.fields['odor_description'].widget = forms.Textarea(attrs=None)
        self.fields['odor_description'].required = True
        self.fields['odor_categories'].required = True
