# Generated by Django 2.0.4 on 2018-10-26 22:15

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0056_bioactive_approval_date'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compoundsource',
            name='product_number',
        ),
        migrations.AlterField(
            model_name='activity',
            name='classification',
            field=models.CharField(choices=[('AT', 'Alimentary tract and metabolism'), ('AI', 'Antiinfectives'), ('AN', 'Anticancer agents'), ('AP', 'Antiparasitics'), ('CV', 'Cardiovascular system'), ('DM', 'Dermatologicals'), ('GU', 'Genito-urinary and sex hormones'), ('MI', 'Miscellaneous'), ('MS', 'Musculo-skeletal system'), ('NS', 'Nervous system'), ('RS', 'Respiratory system')], db_index=True, max_length=2),
        ),
    ]
