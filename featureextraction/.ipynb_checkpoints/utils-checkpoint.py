from . import aminos
import os,sys,inspect
import pandas as pd

parentdir = os.path.dirname(os.getcwd())
sys.path.insert(0,parentdir) 

data = pd.read_csv('cleaned_data.csv')


def get_List_of_Families(data):
    families = data.groupby(['family_name']).agg(
        ['nunique']).reset_index(drop=False)
    families.columns = ['family_name', '#sequences']
    return families


def get_Frequencies(sequence):
    count_freq = {}

    for amino in aminos.amino_acids:
        count_freq[amino] = 0

    for amino in sequence:
        count_freq[amino] += 1

    return count_freq

def getMolecularProperties(sequence):

    pos_MW = 0
    pos_IW = 0
    
    for pos, amino in enumerate(sequence):
        
        #finding Positional Molecular Weight
        pos_MW += molproperties.findPositionalMolecularWeight(amino, pos+1)
        
        #finding Positional Isoelectric Weight
        pos_IW += molproperties.findPositionalIsoelectricWeight(amino, pos+1)
    
    length = len(row['Sequence'])
    
    return (pos_MW / length, pos_IW / length)



def createDictionary(min_n_value, max_n_value):
    extracted_features = {}
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
    
    for amino in amino_acids:
        extracted_features[f'{amino}-mean'] = []
        extracted_features[f'{amino}-SD'] = []
        extracted_features[f'{amino}-CV'] = []
    
    return extracted_features



def extractFeatures(min_n_value = 2, max_n_value = 5):
    extracted_features = createDictionary(min_n_value, max_n_value)
    
    for index, row in data.iterrows():
        if index == 2:
            break
        
        pamw, paiw = utils.getMolecularProperties(row['Sequence'], extracted_features)
        extracted_features['PAMW'].append(pamw)
        extracted_features['PAMW'].append(paiw)
        
        freq_Hpho, freq_Hphi, freq_Neu = hydropathy.findHydropathyComposition(row, index)
        extracted_features['%Hpho'].append(HydropathyComp['Hpho'])
        extracted_features['%Hphi'].append(HydropathyComp['Hphi'])
        extracted_features['%Neu'].append(HydropathyComp['Neu'])
        
        freq_Relative_Hydropathy = hydropathy.findHydropathyTransmission(row['Sequence'])
        extracted_features['%Hpho_Hpho'].append(freq_Relative_Hydropathy['Hpho_Hpho'])
        extracted_features['%Hpho_Hphi'].append(freq_Relative_Hydropathy['Hpho_Hphi'])
        extracted_features['%Hpho_Neu'].append(freq_Relative_Hydropathy['Hpho_Neu'])

        extracted_features['%Hphi_Hpho'].append(freq_Relative_Hydropathy['Hphi_Hpho'])
        extracted_features['%Hphi_Hphi'].append(freq_Relative_Hydropathy['Hphi_Hphi'])
        extracted_features['%Hphi_Neu'].append(freq_Relative_Hydropathy['Hphi_Neu'])

        extracted_features['%Neu_Hpho'].append(freq_Relative_Hydropathy['Neu_Hpho'])
        extracted_features['%Neu_Hphi'].append(freq_Relative_Hydropathy['Neu_Hphi'])
        extracted_features['%Neu_Neu'].append(freq_Relative_Hydropathy['Hphi_Neu'])
        
                                              
        for n in range(min_n_value, max_n_value+1):
            mean, sd, cv = ngrams.find_n_gram(n, row['Sequence'])
            extracted_features[f'{n}-gram-mean'].append(mean)
            extracted_features[f'{n}-gram-sd'].append(sd)
            extracted_features[f'{n}-gram-cv'].append(cv)
        
        
        sigmoidmean, sigmoidsd, sigmoidcv = sixletter.perform6LetterExchange(row['Sequence'])
        extracted_features['6_letter_mean'].append(sigmoidmean)
        extracted_features['6_letter_sd'].append(sigmoidsd)
        extracted_features['6_letter_cv'].append(sigmoidcv)
        
        
        n_mean, n_SD, n_CV = distanceencoding.performDistanceBasedEncoding(sequence)
        for amino in amino_acids:
            extracted_features[f'{amino}-mean'].append(n_mean[f'{amino}'])
            extracted_features[f'{amino}-SD'].append(n_SD[f'{amino}'])
            extracted_features[f'{amino}-CV'].append(n_CV[f'{amino}'])
            
    return extracted_features
    