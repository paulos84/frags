# Generated by Django 2.0.4 on 2018-11-19 23:05

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0071_auto_20181115_2351'),
    ]

    operations = [
        migrations.AddField(
            model_name='substructure',
            name='excludes',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=2), blank=True, default=list, help_text='SMILES excludes element chars', size=None),
        ),
        migrations.AddField(
            model_name='substructure',
            name='includes',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=2), blank=True, default=list, help_text='SMILES must contain element chars', size=None),
        ),
        migrations.AddField(
            model_name='substructure',
            name='max_carbons',
            field=models.IntegerField(blank=True, max_length=2, null=True),
        ),
    ]
