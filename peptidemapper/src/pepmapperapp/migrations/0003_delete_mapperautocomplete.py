# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('pepmapperapp', '0002_mapperautocomplete'),
    ]

    operations = [
        migrations.DeleteModel(
            name='MapperAutoComplete',
        ),
    ]
