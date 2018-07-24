from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import ugettext_lazy as _

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


class UserOdorantSourceCreateForm(forms.ModelForm):

    class Meta:
        model = CompoundSource
        fields = ['price', 'currency', 'amount', 'specification', 'supplier', 'product_number', 'url']
        labels = {
            'url': 'Webpage URL',
        }



#
#
#
# from django import forms
# from django.core.exceptions import ObjectDoesNotExist, ValidationError
#
# from compounds.models import Odorant
#
#
# class OdorantCreateForm(forms.ModelForm):
#
#     """ Form for creating Compound model instances """
#
#     smiles = forms.CharField(widget=forms.HiddenInput())
#     iupac_name = forms.CharField(widget=forms.HiddenInput())
#     cid_number = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'hidden_cid'}))
#
#     class Meta:
#         model = Odorant
#         fields = ['cas_number', 'odor_description', 'odor_categories', 'trade_name', ]
#         widgets = {
#             'odor_description': forms.Textarea(attrs={'rows': 2, 'cols': 42, }),
#             'cas_number': forms.TextInput(attrs={'style': 'border-color: green;', 'size': 44,
#                                                  'placeholder': 'e.g. 80-54-6', }),
#             'odor_categories': forms.SelectMultiple(attrs={'size': '6', }),
#             'trade_name': forms.TextInput(attrs={'size': 40, }),
#         }
#
#     def __init__(self, *args, **kwargs):
#         # self.request = kwargs.pop('request', None)
#         super(OdorantCreateForm, self).__init__(*args, **kwargs)
#         self.fields['odor_description'].required = True
#         self.fields['odor_categories'].required = True
#
#     def clean_cas_number(self):
#         cas_no = self.cleaned_data.get('cas_number')
#         try:
#             Odorant.objects.get(cas_number__exact=cas_no)
#             raise ValidationError('Compound already exists in database')
#         except ObjectDoesNotExist:
#             return cas_no
#
#     def save(self, commit=True):
#         obj = super(OdorantCreateForm, self).save(commit=False)
#         obj.iupac_name = self.cleaned_data['iupac_name']
#         obj.smiles = self.cleaned_data['smiles']
#         obj.cid_number = self.cleaned_data['cid_number']
#         # if self.request:
#         #     obj.created_by = self.request.user
#         if commit:
#             obj.save()
#         return obj
