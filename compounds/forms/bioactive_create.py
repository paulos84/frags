from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError

from compounds.models import Bioactive


class OdorantCreateForm(forms.ModelForm):

    """ Form for creating Compound model instances """

    smiles = forms.CharField(widget=forms.HiddenInput())
    iupac_name = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'iupac_name_field_id'}))
    cid_number = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'hidden_cid'}))

    class Meta:
        model = Bioactive
        fields = ['cas_number', 'odor_description', 'odor_categories', 'trade_name', ]
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 2, 'cols': 42, }),
            'cas_number': forms.TextInput(attrs={'style': 'border-color: green;', 'size': 44,
                                                 'placeholder': 'e.g. 80-54-6',
                                                 'id': 'cas_number'}),
            'odor_categories': forms.SelectMultiple(attrs={'size': '6', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
        }