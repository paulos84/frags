# Generated by Django 2.0.4 on 2018-07-24 21:12

import django.contrib.postgres.fields
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0020_auto_20180724_2054'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserBioactive',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('notes', models.TextField(blank=True, max_length=500, null=True)),
                ('literature_refs', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=100), blank=True, null=True, size=None)),
                ('compound', models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, to='compounds.Bioactive')),
                ('user', models.ForeignKey(blank=True, on_delete=django.db.models.deletion.CASCADE, to='compounds.Profile')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
