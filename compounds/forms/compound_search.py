from django import forms


class OdorantSearchForm(forms.Form):

    cas_number = forms.CharField(
        label='CAS number',
        max_length=14,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. 80-54-6', }),
    )
    chemical_name = forms.CharField(
        label='Chemical name',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. Thujone', }),
    )
    iupac_name = forms.CharField(
        label='IUPAC name contains',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. 2-phenylethyl', }),
    )


class BioactiveSearchForm(forms.Form):

    chemical_name = forms.CharField(
        label='Chemical name',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. aspirin', }),
    )
    iupac_name = forms.CharField(
        label='IUPAC name contains',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. 2-phenylethyl', }),
    )
    inchikey = forms.CharField(
        label='InChIKey',
        max_length=45,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N', }),
    )
