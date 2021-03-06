# Generated by Django 2.0.4 on 2018-07-20 23:08

import django.contrib.postgres.fields
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0010_auto_20180721_0028'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserOdorant',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('notes', models.TextField(blank=True, max_length=500, null=True)),
                ('literature_refs', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=100), blank=True, null=True, size=None)),
                ('odorant', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='odorant_notes_set', to='compounds.Odorant')),
                ('user', models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, related_name='notes_set', to='compounds.Profile')),
            ],
        ),
        migrations.RemoveField(
            model_name='usercompound',
            name='compound',
        ),
        migrations.RemoveField(
            model_name='usercompound',
            name='user',
        ),
        migrations.AlterField(
            model_name='usersource',
            name='compound',
            field=models.ManyToManyField(related_name='sources', to='compounds.UserOdorant', verbose_name='User compounds sources'),
        ),
        migrations.DeleteModel(
            name='UserCompound',
        ),
    ]
