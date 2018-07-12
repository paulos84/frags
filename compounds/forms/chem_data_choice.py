from django import forms

from compounds.utils.general import chemical_properties_label_map


class ChemDataChoiceForm(forms.Form):

    property_choice = forms.ChoiceField(
        choices=tuple(chemical_properties_label_map.items()),
        label='Select chemical property statistics',
        )
