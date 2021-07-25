from . import aminos, distanceencoding, hydropathy, molproperties, ngrams, sixletter, disulphide
import os
import sys
import inspect
import pandas as pd

parentdir = os.path.dirname(os.getcwd())
sys.path.insert(0, parentdir)

data = pd.read_csv('cleaned_data.csv')

amino_acids = aminos.amino_acids
properties = aminos.properties


def getMolecularProperties(sequence):

    pos_MW = 0
    pos_IW = 0

    for pos, amino in enumerate(sequence):

        # finding Positional Molecular Weight
        pos_MW += molproperties.findPositionalMolecularWeight(amino, pos+1)

        # finding Positional Isoelectric Weight
        pos_IW += molproperties.findPositionalIsoelectricWeight(amino, pos+1)

    length = len(sequence)

    return (pos_MW / length, pos_IW / length)


def createDictionary(min_n_value, max_n_value, option):
    extracted_features = {}
    
    if option == 'a' or option == 'b':
        extracted_features['PAMW'] = []
        extracted_features['PAIW'] = []

        extracted_features['%Hpho'] = []
        extracted_features['%Hphi'] = []
        extracted_features['%Neu'] = []

        extracted_features['%Hpho_Hpho'] = []
        extracted_features['%Hpho_Hphi'] = []
        extracted_features['%Hpho_Neu'] = []

        extracted_features['%Hphi_Hpho'] = []
        extracted_features['%Hphi_Hphi'] = []
        extracted_features['%Hphi_Neu'] = []

        extracted_features['%Neu_Hpho'] = []
        extracted_features['%Neu_Hphi'] = []
        extracted_features['%Neu_Neu'] = []

        extracted_features['6_letter_mean'] = []
        extracted_features['6_letter_sd'] = []
        extracted_features['6_letter_cv'] = []

        for n in range(min_n_value, max_n_value+1):
            extracted_features[f'{n}-gram-mean'] = []
            extracted_features[f'{n}-gram-sd'] = []
            extracted_features[f'{n}-gram-cv'] = []
            
            
    if option == 'a' or option == 's':
        for amino in amino_acids:
            extracted_features[f'{amino}-mean'] = []
            extracted_features[f'{amino}-SD'] = []
            extracted_features[f'{amino}-CV'] = []

        extracted_features['para_ds_mean'] = []
        extracted_features['alt_ds_mean'] = []
        extracted_features['quad_ds_mean'] = []

    extracted_features['family_name'] = []

    return extracted_features


def extractFeatures(limit=0, min_n_value=2, max_n_value=5, option = 'a'):
    extracted_features = createDictionary(min_n_value, max_n_value, option)

    for index, row in data.iterrows():

        if limit != 0:
            if index == limit:
                break
        
        if option == 'a' or option == 'b':
            # behavioural features
            pamw, paiw = getMolecularProperties(
                row['Sequence'])
            extracted_features['PAMW'].append(pamw[0])
            extracted_features['PAIW'].append(paiw[0])

            freq_Hpho, freq_Hphi, freq_Neu = hydropathy.findHydropathyComposition(
                row, index)
            extracted_features['%Hpho'].append(freq_Hpho)
            extracted_features['%Hphi'].append(freq_Hphi)
            extracted_features['%Neu'].append(freq_Neu)

            freq_Relative_Hydropathy = hydropathy.findHydropathyTransmission(
                row['Sequence'])
            extracted_features['%Hpho_Hpho'].append(
                freq_Relative_Hydropathy['Hpho_Hpho'])
            extracted_features['%Hpho_Hphi'].append(
                freq_Relative_Hydropathy['Hpho_Hphi'])
            extracted_features['%Hpho_Neu'].append(
                freq_Relative_Hydropathy['Hpho_Neu'])

            extracted_features['%Hphi_Hpho'].append(
                freq_Relative_Hydropathy['Hphi_Hpho'])
            extracted_features['%Hphi_Hphi'].append(
                freq_Relative_Hydropathy['Hphi_Hphi'])
            extracted_features['%Hphi_Neu'].append(
                freq_Relative_Hydropathy['Hphi_Neu'])

            extracted_features['%Neu_Hpho'].append(
                freq_Relative_Hydropathy['Neu_Hpho'])
            extracted_features['%Neu_Hphi'].append(
                freq_Relative_Hydropathy['Neu_Hphi'])
            extracted_features['%Neu_Neu'].append(
                freq_Relative_Hydropathy['Hphi_Neu'])

            for n in range(min_n_value, max_n_value+1):
                mean, sd, cv = ngrams.find_n_gram(n, row['Sequence'])
                extracted_features[f'{n}-gram-mean'].append(mean)
                extracted_features[f'{n}-gram-sd'].append(sd)
                extracted_features[f'{n}-gram-cv'].append(cv)

            sigmoidmean, sigmoidsd, sigmoidcv = sixletter.perform6LetterExchange(
                row['Sequence'])
            extracted_features['6_letter_mean'].append(sigmoidmean)
            extracted_features['6_letter_sd'].append(sigmoidsd)
            extracted_features['6_letter_cv'].append(sigmoidcv)
            
        
        if option == 'a' or option == 's':
            # structural features
            n_mean, n_SD, n_CV = distanceencoding.performDistanceBasedEncoding(
                row['Sequence'])
            for amino in amino_acids:
                extracted_features[f'{amino}-mean'].append(n_mean[f'{amino}'])
                extracted_features[f'{amino}-SD'].append(n_SD[f'{amino}'])
                extracted_features[f'{amino}-CV'].append(n_CV[f'{amino}'])

            para_mean, alt_mean, quad_mean = disulphide.getDisulphideProperties(
                row['Sequence'])
            extracted_features['para_ds_mean'].append(para_mean)
            extracted_features['alt_ds_mean'].append(alt_mean)
            extracted_features['quad_ds_mean'].append(quad_mean)

        extracted_features['family_name'].append(row['family_name'])

    return extracted_features
