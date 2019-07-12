""" Generate BcForms for all of the protein complexes in ComplexPortal

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-07-12
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bcforms
import bpforms
import csv
import os
import re
import requests

IN_FILENAME = os.path.join('examples', 'complextab_homo_sapiens.tsv')
OUT_FILENAME = os.path.join('examples', 'complextab_homo_sapiens_bcforms.tsv')

RNACENTRAL_URL = 'https://rnacentral.org/api/v1/rna/{}/'

def run():
    """ read all human protein complex from complexportal.
        If the subunits are well-defined, represent them in BcForm

    """
    complexes_all = read_from_tsv(IN_FILENAME)
    # print(len(complexes_all))

    complexes_valid = get_valid_complexes(complexes_all)

    for complex in complexes_valid:
        complex['composition'] = get_bcform(complex['composition'])


def read_from_tsv(filename):
    """ Extract complex id and complex composition from TSV file

    Args:
        filename (:obj:`str`): filename of the tsv file

    Returns:
        :obj:`list` of :obj:`dict`: list of complexes
    """
    complexes = []
    with open(filename, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            subunits = row['Identifiers (and stoichiometry) of molecules in complex'].split('|')
            composition = []
            for subunit in subunits:
                words = subunit.split('(')
                composition.append({'id' : words[0], 'stoichiometry':int(words[1][:-1])})
            complex = {'id': row['#Complex ac'], 'composition': composition}
            complexes.append(complex)

    return complexes

def get_valid_complexes(complexes):
    """ get the subset of complexes where all stoichiometry is known

    Args:
        complexes (:obj:`list` of :obj:`dict`): list of all complexes

    Returns:
        :obj:`list` of :obj:`dict`: list of valid complexes
    """
    complexes_valid = []

    for complex in complexes:
        composition = complex['composition']
        valid = True
        for subunit in composition:
            if subunit['stoichiometry'] == 0:
                valid = False
                break

        if valid:
            complexes_valid.append(complex)

    return complexes_valid

def get_bcform(composition):
    """ create BcForm object for the complex

    Args:
        composition (:obj:`list` of :obj:`dict`): list of all subunits the complex consists of

    Returns:
        :obj:`BcForm`: BcForm representation of the complex composition
    """

    bc_form = bcforms.BcForm().from_set(composition)
    # print(str(bc_form))

    # fetch structure information
    for subunit in bc_form.subunits:
        id = subunit.id
        if id.startswith('CHEBI:'):
            # if subunit is small chemical, get structure from CHEBI
            # bc_form.set_subunit_attribute(id, 'structure', )
            pass
        elif re.match(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", id):
            # if subunit is protein, get sequence from UniProt
            # bc_form.set_subunit_attribute(id, 'structure', )
            pass
        elif re.match(r"URS[0-9A-F]{10}", id):
            # if subunit is rna, get sequence from RNAcentral
            response = requests.get(url=RNACENTRAL_URL.format(id))
            if response.status_code != 200:
                return None
            fasta = response.json()['sequence']
            bc_form.set_subunit_attribute(id, 'structure', bpforms.RnaForm().from_str(fasta))

        elif id.startswith('CPX-'):
            # if subunit is a complex, because bcforms do not fully support
            # nested complex, returning None
            # print('BcForm does not fully support nested complexes, returning None as BcForm')
            return None
        else:
            # print('Subunit id not recognized, returning None as BcForm')
            return None

    return bc_form

if __name__ == '__main__':
    run()
