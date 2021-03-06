# Generated by Django 2.0.4 on 2018-07-21 17:06

import django.core.validators
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0014_auto_20180721_1855'),
    ]

    operations = [
        migrations.AlterField(
            model_name='userodorantsource',
            name='currency',
            field=models.CharField(choices=[('USD', 'US Dollars'), ('HKD', 'Hong Kong Dollar'), ('CNY', 'Chinese Yuan'), ('JPY', 'Japanese Yen'), ('GBP', 'British Pound'), ('EUR', 'Euro'), ('IND', 'Indian Rupee')], max_length=3, validators=[django.core.validators.RegexValidator('(?<![A-Z])[A-Z]{3}(?![A-Z])', 'Format must be e.g. USD')]),
        ),
        migrations.AlterField(
            model_name='userodorantsource',
            name='specification',
            field=models.CharField(max_length=100),
        ),
        migrations.AlterField(
            model_name='userodorantsource',
            name='supplier',
            field=models.CharField(max_length=50),
        ),
    ]
