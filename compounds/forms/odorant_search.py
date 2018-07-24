from django import forms


class OdorantSearchForm(forms.Form):

    cas_number = forms.CharField(
        label='CAS number',
        max_length=12,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. 80-54-6', }),
    )
    iupac_name = forms.CharField(
        label='Name contains',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. 2-phenylethyl', }),
    )


class BioactiveSearchForm(forms.Form):

    inchikey = forms.CharField(
        label='InChIKey',
        max_length=45,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N', }),
    )
    iupac_name = forms.CharField(
        label='Name contains',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. 2-phenylethyl', }),
    )