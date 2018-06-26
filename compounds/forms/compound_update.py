from django import forms

from compounds.models import Compound


class CompoundUpdateForm(forms.ModelForm):

    """ Form for updating Compound model instances """

    class Meta:
        model = Compound
        fields = ('odor_description', 'odor_categories', 'trade_name')
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 2, 'cols': 42, }),
            'odor_categories': forms.SelectMultiple(attrs={'size': '6', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
        }
