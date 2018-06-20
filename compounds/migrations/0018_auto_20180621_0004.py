# Generated by Django 2.0.4 on 2018-06-20 22:04

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0017_auto_20180618_2119'),
    ]

    operations = [
        migrations.AddField(
            model_name='substructure',
            name='slug',
            field=models.SlugField(default=''),
        ),
        migrations.AlterField(
            model_name='compoundnotes',
            name='compound',
            field=models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, related_name='notes_set', to='compounds.Compound'),
        ),
        migrations.AlterField(
            model_name='compoundnotes',
            name='user',
            field=models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, related_name='notes_set', to='compounds.Profile'),
        ),
        migrations.AlterField(
            model_name='substructure',
            name='name',
            field=models.CharField(max_length=50, verbose_name='Substructure class name'),
        ),
    ]