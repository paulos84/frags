from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

from compounds.models import UserBioactive


class CompoundNotesForm(forms.Form):
    """
    Form for users to create or update a CompoundNote instance
    """
    notes = forms.CharField(
        max_length=200,
        widget=forms.Textarea(
            attrs={'rows': 3, 'cols': 30, 'placeholder': 'Enter notes',
                   'style': 'border-color: green;', }),
    )

    def __init__(self, *args, **kwargs):
        user_auth = kwargs.pop('user_auth', None)
        notes = kwargs.pop('notes', None)
        super(CompoundNotesForm, self).__init__(*args, **kwargs)
        if notes:
            self.fields['notes'].initial = notes
        if not user_auth:
            self.fields['notes'].widget.attrs['placeholder'] = 'Login to access notes'


class UserBioactiveChemDataForm(forms.Form):
    """
    Form for users to update chemical_data field on a UserBioactive instance
    """
    user_bioactive = forms.ModelChoiceField(
        queryset=UserBioactive.objects.all(),
        widget=forms.HiddenInput(),
        required=False
    )
    label = forms.CharField(
        max_length=25,
        label='',
        required=True,
        widget=forms.TextInput(attrs={'placeholder': 'Label', })
    )
    value = forms.CharField(
        max_length=30,
        label='',
        required=True,
        widget=forms.TextInput(attrs={'placeholder': 'Data', })
    )

    def save(self):
        instance = self.cleaned_data['user_bioactive']
        json_data = {self.cleaned_data['label']: self.cleaned_data['value']}
        instance.chemical_data.update(json_data)
        instance.save()


class UserLiteratureRefsForm(forms.Form):
    lit_ref_numbers = forms.MultipleChoiceField(
        widget=forms.CheckboxSelectMultiple,
        required=True,
    )

    def __init__(self, *args, **kwargs):
        available_choices = kwargs.pop('lit_records')
        super(UserLiteratureRefsForm, self).__init__(*args, **kwargs)
        self.fields['lit_ref_numbers'].choices = [(a, '') for a in available_choices]


class SignupForm(UserCreationForm):
    email = forms.EmailField(max_length=200, help_text='Required')

    class Meta:
        model = User
        fields = ('username', 'email', 'password1', 'password2')

    def __init__(self, *args, **kwargs):
        super(UserCreationForm, self).__init__(*args, **kwargs)
        for field in ('username', 'email', 'password1', 'password2'):
            self.fields[field].widget.attrs['class'] = 'form-control'
            self.fields[field].help_text = ''


class ContactForm(forms.Form):
    contact_email = forms.EmailField(required=True)
    subject = forms.CharField(required=True)
    message = forms.CharField(widget=forms.Textarea, required=True)
