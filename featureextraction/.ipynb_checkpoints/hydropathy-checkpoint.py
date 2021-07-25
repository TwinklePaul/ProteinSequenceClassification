from . import aminos

properties = aminos.properties


def findHydropathyComposition(data_row, index):
    freq_Hpho = 0
    freq_Hphi = 0
    freq_Neu = 0

    for amino in aminos.amino_acids:
        freq_amino_in_sequence = data_row[amino]
        amino_hydropathy_type = properties.loc[properties['Tag']
                                               == amino]['Hydropathy_Label'].values

        if amino_hydropathy_type < 0:
            freq_Hpho += freq_amino_in_sequence
        elif amino_hydropathy_type > 0:
            freq_Hphi += freq_amino_in_sequence
        else:
            freq_Neu += freq_amino_in_sequence

    freq_Hpho = (freq_Hpho * 100) / length
    freq_Hphi = (freq_Hphi * 100) / length
    freq_Neu = (freq_Neu * 100) / length
    
    return (freq_Hpho, freq_Hphi, freq_Neu)


def findHydropathyTransmission(sequence):

    freq_Relative_Hydropathy = {
        'Hpho_Hpho': 0,
        'Hpho_Hphi': 0,
        'Hpho_Neu': 0,
        'Hphi_Hpho': 0,
        'Hphi_Hphi': 0,
        'Hphi_Neu': 0,
        'Neu_Hpho': 0,
        'Neu_Hphi': 0,
        'Neu_Neu': 0
    }

    for pos, amino in enumerate(sequence):

        if pos == len(sequence) - 1:
            break

        next_amino = sequence[pos+1]

        amino_hydropathy_type = properties.loc[properties['Tag']
                                               == amino]['Hydropathy_Label'].values
        next_amino_hydropathy_type = properties.loc[properties['Tag']
                                                    == next_amino]['Hydropathy_Label'].values

        if amino_hydropathy_type < 0:

            if next_amino_hydropathy_type < 0:
                freq_Relative_Hydropathy['Hpho_Hpho'] += 1
            elif next_amino_hydropathy_type > 0:
                freq_Relative_Hydropathy['Hpho_Hphi'] += 1
            else:
                freq_Relative_Hydropathy['Hpho_Neu'] += 1

        if amino_hydropathy_type > 0:

            if next_amino_hydropathy_type < 0:
                freq_Relative_Hydropathy['Hphi_Hpho'] += 1
            elif next_amino_hydropathy_type > 0:
                freq_Relative_Hydropathy['Hphi_Hphi'] += 1
            else:
                freq_Relative_Hydropathy['Hphi_Neu'] += 1

        if amino_hydropathy_type == 0:

            if next_amino_hydropathy_type < 0:
                freq_Relative_Hydropathy['Neu_Hpho'] += 1
            elif next_amino_hydropathy_type > 0:
                freq_Relative_Hydropathy['Neu_Hphi'] += 1
            else:
                freq_Relative_Hydropathy['Neu_Neu'] += 1
                
    freq_Relative_Hydropathy['Hpho_Hpho'] = (freq_Relative_Hydropathy['Hpho_Hpho'] * 100) / length
    freq_Relative_Hydropathy['Hpho_Hphi'] = (freq_Relative_Hydropathy['Hpho_Hphi'] * 100) / length
    freq_Relative_Hydropathy['Hpho_Neu'] = (freq_Relative_Hydropathy['Hpho_Neu'] * 100) / length
    
    freq_Relative_Hydropathy['Hphi_Hpho'] = (freq_Relative_Hydropathy['Hphi_Hpho'] * 100) / length
    freq_Relative_Hydropathy['Hphi_Hphi'] = (freq_Relative_Hydropathy['Hphi_Hphi'] * 100) / length
    freq_Relative_Hydropathy['Hphi_Neu'] = (freq_Relative_Hydropathy['Hphi_Neu'] * 100) / length
    
    freq_Relative_Hydropathy['Neu_Hpho'] = (freq_Relative_Hydropathy['Neu_Hpho'] * 100) / length
    freq_Relative_Hydropathy['Neu_Hphi'] = (freq_Relative_Hydropathy['Neu_Hphi'] * 100) / length
    freq_Relative_Hydropathy['Neu_Neu'] = (freq_Relative_Hydropathy['Neu_Neu'] * 100) / length


    return freq_Relative_Hydropathy
