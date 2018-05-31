from django import forms
from django.core.exceptions import ObjectDoesNotExist, ValidationError
from django.db.models import Q

from compounds.models import Compound


class CompoundCreateForm(forms.ModelForm):
    smiles = forms.CharField(widget=forms.HiddenInput())
    iupac_name = forms.CharField(widget=forms.HiddenInput())
    cid_number = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'hidden_cid'}))

    class Meta:
        model = Compound
        fields = ['cas_number', 'odor_description', 'odor_category', 'trade_name', 'supplier']
        widgets = {
            'odor_description': forms.Textarea(attrs={'rows': 5, 'cols': 52, }),
            'cas_number': forms.TextInput(attrs={'style': 'border-color: green;', 'size': 50,
                                                 'placeholder': 'Enter CAS number to lookup compound', }),
            'odor_category': forms.SelectMultiple(attrs={'size': '8', }),
            'trade_name': forms.TextInput(attrs={'size': 40, }),
        }

    def __init__(self, *args, **kwargs):
        super(CompoundCreateForm, self).__init__(*args, **kwargs)
        self.fields['odor_description'].required = True
        self.fields['odor_category'].required = True

    def clean_cas_number(self):
        cas_no = self.cleaned_data.get('cas_number')
        try:
            Compound.objects.get(
                Q(cas_number__exact=cas_no) | Q(additional_cas__contains=cas_no)
            )
            raise ValidationError('Compound already exists in database')
        except ObjectDoesNotExist:
            return cas_no

    def save(self, commit=True):
        obj = super(CompoundCreateForm, self).save(commit=False)
        obj.iupac_name = self.cleaned_data['iupac_name']
        obj.smiles = self.cleaned_data['smiles']
        obj.cid_number = self.cleaned_data['cid_number']
        if commit:
            obj.save()
        return obj


class CompoundUpdateForm(forms.ModelForm):

    # https: // docs.djangoproject.com / en / 2.0 / topics / forms / modelforms /  # providing-initial-values


    class Meta:
        model = Compound
        fields = ['odor_description', 'odor_category', 'supplier', 'trade_name']

    def __init__(self, *args, **kwargs):
        super(CompoundUpdateForm, self).__init__(*args, **kwargs)
        self.fields['odor_description'].widget = forms.Textarea(attrs=None)
        self.fields['odor_description'].required = True
        self.fields['odor_category'].required = True
        # self.fields['name'].widget.attrs.update({'class': 'special'})
        # self.fields['comment'].widget.attrs.update(size='40')



#sharegp = forms.ModelChoiceField(label='Share with groups', queryset=usergroups)

    # def clean_zip_code(self):
    #     return self.cleaned_data['zip_code'].upper()

    # def clean_account_type(self):
    #     account_type = self.cleaned_data["account_type"]
    #     if account_type == "select":
    #         raise forms.ValidationError("Select account type.")
    #     return account_type
    """

    def clean(self):
        cleaned_data = super().clean()

        cc_myself = cleaned_data.get("cc_myself")
        subject = cleaned_data.get("subject")

        if cc_myself and subject and "help" not in subject:
            msg = "Must put 'help' in subject when cc'ing yourself."
            self.add_error('cc_myself', msg)
            self.add_error('subject', msg)

Field for form use a difference library to create a from. You need to import django.forms and use form.XXX for specific Field

from django import forms


class StudentForm(ModelForm):
    class Meta:
        model = Student

    subject = forms.CharField(label='New label')

In order to customize field in model form, you don't need to create it manually. Django model fields have special attributes:

    verbose_name (goes to label of the field)
    help_text (by default rendered as additional description below the field)

So, all you need is:

class Student(models.Model):
    name = models.CharField(max_length=40,
                            verbose_name="Student's Name",
                            help_text="Please tell me your name")  # Optional
    last_name = models.CharFIeld(max_length=40)
    ...

Then you don't need to do any customization in model form.
    """