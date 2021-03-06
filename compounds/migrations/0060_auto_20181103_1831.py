# Generated by Django 2.0.4 on 2018-11-03 18:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0059_odorant_hit_count'),
    ]

    operations = [
        migrations.AddField(
            model_name='bioactive',
            name='notes',
            field=models.CharField(blank=True, max_length=500, null=True, verbose_name='extra information'),
        ),
        migrations.AlterField(
            model_name='activity',
            name='classification',
            field=models.CharField(blank=True, choices=[('AT', 'Alimentary tract and metabolism'), ('AI', 'Antiinfectives'), ('AN', 'Anticancer agents'), ('AP', 'Antiparasitics'), ('CV', 'Cardiovascular system'), ('DM', 'Dermatologicals'), ('GU', 'Genito-urinary and sex hormones'), ('MI', 'Miscellaneous'), ('MS', 'Musculo-skeletal system'), ('NS', 'Nervous system'), ('RS', 'Respiratory system')], db_index=True, max_length=2, null=True),
        ),
    ]
