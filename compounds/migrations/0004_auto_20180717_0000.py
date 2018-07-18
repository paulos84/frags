# Generated by Django 2.0.4 on 2018-07-16 22:00

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0003_auto_20180717_0000'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usernotes',
            name='notes',
            field=models.TextField(blank=True, max_length=500, null=True),
        ),
        migrations.AlterField(
            model_name='usernotes',
            name='user',
            field=models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, related_name='notes_set', to='compounds.Profile'),
        ),
    ]