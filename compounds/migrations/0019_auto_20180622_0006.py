# Generated by Django 2.0.4 on 2018-06-21 22:06

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0018_auto_20180621_0004'),
    ]

    operations = [
        migrations.AlterField(
            model_name='substructure',
            name='description',
            field=models.CharField(blank=True, default='', max_length=1000),
        ),
        migrations.AlterField(
            model_name='substructure',
            name='slug',
            field=models.SlugField(blank=True, default=''),
        ),
    ]
