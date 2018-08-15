from django import forms
from django.core.exceptions import ValidationError

from compounds.models import CompoundSource


def validate_file_extension(value):
    if not value.name.endswith('.csv'):
        raise ValidationError('Must be a csv file')


class UserSourceCsvUploadForm(forms.Form):
    csv_file = forms.FileField(
        validators=[validate_file_extension],
        help_text='Headers ordered as above, no currency',
        required=False,
    )
    currency = forms.ChoiceField(
        choices=CompoundSource.currency_choices,
    )

    def __init__(self, *args, **kwargs):
        super(UserSourceCsvUploadForm, self).__init__(*args, **kwargs)
        self.fields['currency'].widget.attrs['style'] = "width:200px"


class CompoundSourceCreateForm(forms.ModelForm):

    class Meta:
        model = CompoundSource
        fields = ['price', 'currency', 'amount', 'unit', 'specification', 'supplier', 'product_number', 'url']
        labels = {
            'url': 'Webpage URL',
            'unit': 'Units'
        }

    def __init__(self, *args, **kwargs):
        super(CompoundSourceCreateForm, self).__init__(*args, **kwargs)
        for field in self.fields:
            self.fields[field].widget.attrs['style'] = "width:300px"
