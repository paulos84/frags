# Generated by Django 2.0.4 on 2018-11-30 21:01

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0080_auto_20181130_2040'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='development',
            name='lit_refs',
        ),
    ]
