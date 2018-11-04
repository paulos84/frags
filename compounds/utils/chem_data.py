
chemical_properties_label_map = {'xlogp': 'Partition coefficient (xlogp)',
                                 'hac': 'Non-hydrogen atom count',
                                 'rbc': 'Rotable bond count',
                                 'hetac': 'Heteroatom count',
                                 'mw': 'Molecular weight'}

colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd',
          '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
          '#17becf', '#9edae5', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c', '#fdae6b',
          '#fdd0a2', '#31a354', '#74c476', '#a1d99b', '#c7e9c0', '#756bb1', '#9e9ac8', '#bcbddc', '#dadaeb',
          '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']


def dict_from_query_object(smiles, pcp_data, additional):
    chem_dict = {a: getattr(pcp_data, b) for a, b in
                 (('xlogp', 'xlogp'), ('hac', 'heavy_atom_count'), ('rbc', 'rotatable_bond_count'))}
    chem_dict.update({
        'mw': int(pcp_data.molecular_weight),
        'synonyms': ', '.join(pcp_data.synonyms[:5]),
        'hetac': len([i for i in smiles.split('.')[0] if i in
                      ['N', 'S', 'O', 'Cl', 'F', 'Br', 'B', 'P']])
    })
    if additional:
        extra_available = ['h_bond_acceptor_count', 'h_bond_donor_count', 'complexity', 'atom_stereo_count',
                           'bond_stereo_count']
        chem_dict.update(
            {k if k != 'rotable_bond_count' else 'rbc': int(getattr(pcp_data, k))
             for k in extra_available}
        )
    return chem_dict
