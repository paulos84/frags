# Generated by Django 2.0.4 on 2018-06-26 18:43

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0020_auto_20180624_1833'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='CompoundNotes',
            new_name='UserNotes',
        ),
        migrations.AlterField(
            model_name='compound',
            name='created_by',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='compounds', to='compounds.Profile'),
        ),
    ]
