# Generated by Django 2.0.4 on 2018-06-16 23:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0012_substructure'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='substructure',
            name='image',
        ),
        migrations.AddField(
            model_name='substructure',
            name='cid_number',
            field=models.IntegerField(blank=True, default=32, editable=False, verbose_name='PubChem API CID number'),
            preserve_default=False,
        ),
    ]