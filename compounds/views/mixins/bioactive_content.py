from compounds.models import Activity, Bioactive, BioactiveCore


class BioactiveContentMixin:
    """
    Provides logic to enable various bioactive list views to show relevant title and links
    """
    def get_context_data(self, **kwargs):
        context = super(BioactiveContentMixin, self).get_context_data(**kwargs)
        if self.category == 1:
            self.add_pharma_context(context)
        elif self.category == 2:
            self.add_food_context(context)
        label = Bioactive.cat_choices[self.category - 1][1]
        context.update({
            'page_header': label + 's' if not label.endswith('s') else label,
            'category': self.category,
        })
        return context

    def add_pharma_context(self, context):
        context.update({
            'body_systems': Activity.classified_actions_mechs(),
            'drug_actions': Activity.objects.actions().order_by('name'),
        })

    def add_food_context(self, context):
        context.update({
            'substructures': BioactiveCore.objects.food().exclude(name='Oligosaccharides').values('name', 'slug')
        })
