# Generated by Django 2.0.4 on 2018-08-16 22:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0046_auto_20180815_2352'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compoundsource',
            name='amount',
            field=models.FloatField(blank=True, max_length=6, null=True),
        ),
        migrations.AlterField(
            model_name='compoundsource',
            name='price',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=10, null=True),
        ),
        migrations.AlterField(
            model_name='compoundsource',
            name='specification',
            field=models.CharField(blank=True, max_length=30, null=True),
        ),
        migrations.AlterField(
            model_name='compoundsource',
            name='supplier',
            field=models.CharField(blank=True, max_length=25, null=True),
        ),
    ]
