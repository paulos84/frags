from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError

from compounds.models import Bioactive


class BioactiveCreateForm(forms.ModelForm):
    """ Form for creating Compound model instances """

    cas_number = forms.CharField(
        required=False, label='CAS number',
        widget=forms.TextInput(attrs={'size': 44,
                                      'id': 'cas_number_field_id',
                                      'placeholder': 'Retrieve InChIKey (optional)', }
                               ))
    smiles = forms.CharField(
        widget=forms.HiddenInput()
    )
    iupac_name = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'iupac_name_field_id'})
    )
    chemical_name = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'chemical_name_field_id', 'size': 44, })
    )
    cid_number = forms.CharField(
        widget=forms.HiddenInput(attrs={'id': 'hidden_cid'})
    )

    class Meta:
        model = Bioactive
        fields = ['cas_number', 'inchikey', 'category', ]
        widgets = {
            'inchikey': forms.TextInput(attrs={'size': 44, 'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N',
                                               'style': 'border-color: green;', 'id': 'inchikey_field_id'}),
        }
