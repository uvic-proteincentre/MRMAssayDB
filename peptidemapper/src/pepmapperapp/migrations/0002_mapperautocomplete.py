# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('pepmapperapp', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='MapperAutoComplete',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('uniprotacc', models.CharField(max_length=20)),
                ('prot_name', models.CharField(max_length=300)),
                ('gene', models.CharField(max_length=100)),
                ('organism', models.CharField(max_length=100)),
                ('pepseq', models.CharField(max_length=20)),
                ('path_name', models.CharField(max_length=1000)),
                ('dis_mut', models.CharField(max_length=1000)),
                ('go_id', models.CharField(max_length=100)),
                ('go_name', models.CharField(max_length=1000)),
                ('go_term', models.CharField(max_length=100)),
            ],
        ),
    ]
