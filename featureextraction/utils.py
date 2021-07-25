from . import aminos


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
