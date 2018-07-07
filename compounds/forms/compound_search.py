from django import forms


class CompoundSearchForm(forms.Form):

    cas_number = forms.CharField(
        label='CAS number',
        max_length=12,
        required=False,
    )
    iupac_name = forms.CharField(
        label='Compound name',
        max_length=12,
        required=False,
    )
