from django import forms
from django.core.validators import MinLengthValidator

from compounds.models.compound import Compound


class CompoundCreateForm(forms.ModelForm):
    # In instantiated Django forms, fields are kept in a dict-like object. Which means, instead of writing forms in a
    #  way that duplicates the model, a better way is to explicitly modify only what we want to modify:

    class Meta:
        model = Compound
        fields = ['cas_number', 'odor_category', 'odor_description', 'trade_name', 'supplier']

    def __init__(self, *args, **kwargs):
        super(CompoundCreateForm, self).__init__(*args, **kwargs)
        self.fields['cas_number'].help_text = 'Enter a CAS number to retrieve compound',
        self.fields['odor_description'].widget = forms.Textarea(attrs=None)
        self.fields["odor_description"].min_length = 20
        self.fields["odor_description"].validators.append(MinLengthValidator)
        self.fields['odor_description'].required = True
        self.fields['odor_category'].required = True
        # self.fields['odor_description'].widget.attrs.update({'class': 'special'})
        # self.fields['name'].widget.attrs.update({'class': 'special'})
        # self.fields['comment'].widget.attrs.update(size='40')

# TODO: Able to search by drawing in structure also: e.g. https://www.sigmaaldrich.com/catalog/search/substructure/OldSubstructureSearchPage

    # Validation: make sure cas_number not in alternative_cas for all other compounds (i.e. not unique)
    # check that clean applies max_length etc constraints defined in model (unit test)

    def clean(self):
        cleaned_data = super(CompoundCreateForm, self).clean()
        name = cleaned_data.get('name')
        email = cleaned_data.get('email')
        message = cleaned_data.get('message')
        if not name and not email and not message:
            raise forms.ValidationError('You have to write something!')

class CompoundUpdateForm(forms.ModelForm):
    # In instantiated Django forms, fields are kept in a dict-like object. Which means, instead of writing forms in a
    #  way that duplicates the model, a better way is to explicitly modify only what we want to modify:

    class Meta:
        model = Compound
        fields = ['odor_description', 'odor_category', 'supplier', 'trade_name']

    def __init__(self, *args, **kwargs):
        super(CompoundUpdateForm, self).__init__(*args, **kwargs)
        self.fields['odor_description'].widget = forms.Textarea(attrs=None)
        self.fields["odor_description"].min_length = 20
        self.fields["odor_description"].validators.append(MinLengthValidator)
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