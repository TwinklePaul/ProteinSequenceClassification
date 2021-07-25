from . import aminos

properties = aminos.properties


def findPositionalMolecularWeight(amino, pos):
    return properties.loc[properties['Tag'] == amino]['Molecular_Weight'].values * pos


def findPositionalIsoelectricWeight(amino, pos):
    return properties.loc[properties['Tag'] == amino]['IsoElectric_Point'].values * pos
