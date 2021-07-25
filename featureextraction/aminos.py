import pandas as pd

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
               'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


properties = {
    'Tag': amino_acids,

    'Protein_Name': ['Alanine', 'Cysteine', 'Aspartic Acid', 'Glutamic Acid',
                     'Phenylalanine', 'Glycine', 'Histidine', 'Isoleucine',
                     'Lysine', 'Leucine', 'Methionine', 'Asparagine',
                     'Proline', 'Glutamine', 'Arginine', 'Serine',
                     'Threonine', 'Valine', 'Tryptophan', 'Tyrosine'],

    'Molecular_Weight': [89.1, 121.16, 133.11, 147.13, 165.19,
                         75.07, 155.16, 131.18, 146.19, 131.18,
                         149.21, 132.12, 115.13, 146.15, 174.2,
                         105.09, 119.12, 117.15, 204.23, 181.19],

    'IsoElectric_Point': [6, 5.07, 2.77, 3.22, 5.48,
                          5.97, 7.59, 6.02, 9.74, 5.98,
                          5.74, 5.41, 6.3, 5.65, 10.76,
                          5.68, 5.6, 5.96, 5.89, 5.66],

    'Hydropathy_Property': ['Hydrophobic', 'Hydrophobic', 'Hydrophilic', 'Neutral',
                            'Very Hydrophobic', 'Neutral', 'Hydrophilic', 'Very Hydrophobic',
                            'Hydrophilic', 'Very Hydrophobic', 'Very Hydrophobic', 'Neutral',
                            'Hydrophilic', 'Neutral', 'Hydrophilic', 'Neutral',
                            'Neutral', 'Very Hydrophobic', 'Very Hydrophobic', 'Hydrophobic'],

    'Hydropathy_Label': [-1, -1, 1, 0, -2,
                         0, 1, -2, 1, -2,
                         -2, 0, 1, 0, 1,
                         0, 0, -2, -2, -1],

    '6_Letter_Encoding': ['e4', 'e3', 'e2', 'e2', 'e6',
                          'e4', 'e1', 'e5', 'e1', 'e5',
                          'e5', 'e2', 'e4', 'e2', 'e1',
                          'e4', 'e4', 'e5', 'e6', 'e6']
}


properties = pd.DataFrame(properties)
