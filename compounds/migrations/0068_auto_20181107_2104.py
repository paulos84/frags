# Generated by Django 2.0.4 on 2018-11-07 21:04

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0067_enzyme_mechanisms'),
    ]

    operations = [
        migrations.RenameField(
            model_name='enzyme',
            old_name='mechanisms',
            new_name='mechanism',
        ),
        migrations.AlterField(
            model_name='enzyme',
            name='category',
            field=models.IntegerField(choices=[(1, 'Galactosidase'), (2, 'Drug target')], db_index=True, default=2),
        ),
    ]
