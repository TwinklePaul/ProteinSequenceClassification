import math
from . import aminos
from . import ngrams


properties = aminos.properties
count = 0

def perform6LetterEncoding(sequence):
    for amino in sequence:
        exchange_code = properties.loc[properties['Tag'] == amino]['6_Letter_Encoding'].values[0]
        
        sequence = sequence.replace(amino, exchange_code)
    return sequence


def perform6LetterExchange(sequence):

    sequence = perform6LetterEncoding(sequence)

    mean, sd, cv = ngrams.find_n_gram(4, sequence)

    sigmoidmean = 1 / (1 + math.exp(-mean))
    sigmoidsd = 1 / (1 + math.exp(-sd))
    sigmoidcv = 1 / (1 + math.exp(-cv))

    return (sigmoidmean, sigmoidsd, sigmoidcv)
