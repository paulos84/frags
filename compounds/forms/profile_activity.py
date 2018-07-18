from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.contrib.postgres.forms import SimpleArrayField

from compounds.models import UserCompound, Odorant, Profile, UserSource


class CompoundNotesForm(forms.Form):

    """ Form for users to create or update a CompoundNote instance """

    notes = forms.CharField(
        widget=forms.Textarea(
            attrs={'rows': 5, 'cols': 42, 'placeholder': 'Enter notes',
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


class UserLiteratureRefsForm(forms.Form):
    lit_ref_numbers = forms.MultipleChoiceField(
        widget=forms.CheckboxSelectMultiple,
        required=True,
    )

    def __init__(self, *args, **kwargs):
        available_choices = kwargs.pop('lit_records')
        super(UserLiteratureRefsForm, self).__init__(*args, **kwargs)
        self.fields['lit_ref_numbers'].choices = [(a, '') for a in available_choices]


class UserSourcesForm(forms.ModelForm):
    compound = forms.ModelChoiceField(
        queryset=Odorant.objects.all(),
        widget=forms.HiddenInput(),
    )

    class Meta:
        model = UserSource
        fields = ['webpage', 'price_info', 'compound']

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
