from django.db import models
from django.shortcuts import reverse
from django.utils.text import slugify
from django.utils import timezone

from compounds.utils.pipelines_update import AstraZenecaPipeline


class CompanyPipeline(models.Model):

    name = models.CharField(
        max_length=20,
        unique=True
    )
    code = models.CharField(
        max_length=10,
        unique=True,
    )
    info_url = models.URLField(
        max_length=100,
        unique=True,
    )
    updated_at = models.DateField(
        auto_now_add=True,
    )

    def __str__(self):
        return self.name

    @classmethod
    def update_astrazeneca(cls):
        """
        Call using management command, example:
        $ python manage.py update_pipeline --company astrazeneca --settings=frags.settings_LOCAL_10063
        """
        az = AstraZenecaPipeline()
        bioactives = cls.objects.get(name='AstraZeneca').bioactives.all()
        dev_p1 = [b for b in bioactives if b.phase == 1]
        dev_p2 = [b for b in bioactives if b.phase == 2]
        dev_p3 = [b for b in bioactives if b.phase == 3]
        p1_cpds = [a.lower() for a in az.phases['p1']]
        p2_cpds = [a.lower() for a in az.phases['p2']]
        p3_cpds = [a.lower() for a in az.phases['p3']]
        for a in dev_p1:
            if a.bioactive.chemical_name and a.bioactive.chemical_name.lower() not in p1_cpds:
                if a.bioactive.chemical_name.lower() in p2_cpds:
                    a.status = 2
                    a.save()
                    print('Compound {} id:{} moved from phase 1 to phase 2'.format(
                        a.bioactive.chemical_name, a.id))
                else:
                    print('Phase 1 compound {} id:{} not found in p1 or p2 lists'.format(
                        a.bioactive.chemical_name, a.id))
        for a in dev_p2:
            if a.bioactive.chemical_name.lower() not in p2_cpds:
                if a.bioactive.chemical_name.lower() in p3_cpds:
                    a.status = 3
                    a.save()
                    print('Compound {} id:{} moved from phase 2 to phase 3'.format(
                        a.bioactive.chemical_name, a.id))
                else:
                    print('Phase 2 compound {} id:{} not found in p2 or p3 lists'.format(
                        a.bioactive.chemical_name, a.id))
        for a in dev_p3:
            if a.bioactive.chemical_name.lower() not in p3_cpds:
                print('Phase 3 compound {} id:{} not found in p3 list'.format(
                    a.bioactive.chemical_name, a.id))


