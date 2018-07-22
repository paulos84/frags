from django.core.exceptions import ValidationError


chemical_properties_label_map = {'xlogp': 'Partition coefficient(xlogp)',
                                 'hac': 'Non-hydrogen atom count',
                                 'rbc': 'Rotable bond count',
                                 'hetac': 'Heteroatom (N, S, O) count',
                                 'mw': 'Molecular weight'}
