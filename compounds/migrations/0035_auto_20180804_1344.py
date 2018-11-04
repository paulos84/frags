# Generated by Django 2.0.4 on 2018-08-04 11:44

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0034_bioactivecore_related_smiles'),
    ]

    operations = [
        migrations.CreateModel(
            name='Classification',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('category', models.IntegerField(choices=[(1, 'Alimentary tract and metabolism'), (2, 'Antineoplastic and immunomodulating agents'), (3, 'Antiparasitics'), (4, 'Blood modifiers'), (5, 'Cardiovascular system'), (6, 'Dermatologicals'), (7, 'Genito-urinary system and sex hormones'), (8, 'Musculo-skeletal system'), (9, 'Nervous system'), (10, 'Respiratory system'), (11, 'Systemic antiinfectives'), (12, 'Systemic hormones'), (13, 'Various')])),
                ('name', models.CharField(max_length=20, unique=True)),
            ],
        ),
        migrations.AlterField(
            model_name='bioactivecore',
            name='related_smiles',
            field=models.CharField(blank=True, db_index=True, default='', help_text='For substructure matches of close analogs', max_length=100, verbose_name='Close analog SMILES string'),
        ),
        migrations.AddField(
            model_name='bioactive',
            name='classification',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='bioactives', to='compounds.Classification', verbose_name='Primary pharmacological classification'),
        ),
    ]
