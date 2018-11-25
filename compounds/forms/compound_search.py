from django import forms


class ChemNameSearchForm(forms.Form):
    chemical_name = forms.CharField(
        label='Chemical name',
        max_length=40,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. clavulanic acid', }),
    )
    protein_term = forms.CharField(
        label='Drug target protein',
        max_length=40,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'protein or ligand name contains'}),
    )


class ProteinSearchForm(forms.Form):
    search_term = forms.CharField(
        label='Search target proteins',
        max_length=40,
        required=True,
        widget=forms.TextInput(
            attrs={'placeholder': 'protein or ligand name contains'}),
    )


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


class BioactiveSearchForm(forms.Form):

    chemical_name = forms.CharField(
        label='Chemical name',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. clavulanic acid', }),
    )
    iupac_name = forms.CharField(
        label='IUPAC name contains',
        max_length=30,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. azabicyclo', }),
    )
    inchikey = forms.CharField(
        label='InChIKey',
        max_length=45,
        required=False,
        widget=forms.TextInput(
            attrs={'placeholder': 'e.g. YPBKTZBXSBLTDK-PKNBQFBNSA-N', }),
    )
