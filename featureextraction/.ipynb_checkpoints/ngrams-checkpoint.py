import numpy as np


def generatePossiblePatterns(sequence, n):
    list_possible_patterns = []

    for pos in range(0, len(sequence) - n):

        for pattern_length in range(n):
            pattern = sequence[pos: pos+pattern_length+1]

            if pattern not in list_possible_patterns:
                list_possible_patterns.append(pattern)

    return list_possible_patterns


def findOccuranceProbabilities(sequence, n):
    list_occurance_prob = []

    list_possible_patterns = generatePossiblePatterns(sequence, n)

    for pattern in list_possible_patterns:
        match_count = sum(1 for i in range(len(sequence))
                          if sequence.startswith(pattern, i))

        occurance_prob = match_count / (len(sequence) - n - 1)
        list_occurance_prob.append(occurance_prob)

    return list_occurance_prob


def find_n_gram(n, sequence):

    list_occurance_prob = findOccuranceProbabilities(sequence, n)

    mean = np.mean(list_occurance_prob)
    sd = np.std(list_occurance_prob)
    cv = sd/mean

    return (mean, sd, cv)
