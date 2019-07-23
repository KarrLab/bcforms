""" Test of bcforms.rest

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-7-4
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bcforms
from bcforms import core
from bcforms import rest
from wc_utils.util.chem import EmpiricalFormula
import unittest

class RestTestCase(unittest.TestCase):

    def test_PrefixMiddleware(self):
        rest.PrefixMiddleware(rest.app).__call__({'PATH_INFO': 'x'}, lambda x, y: None)


    def test_get_bcform_properties(self):
        client = rest.app.test_client()

        # test validate
        rv = client.post('/api/bcform/', json=dict({
            "form": "abc_a + abc_b"
        }))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b'
        })

        # test when all structure is known (protein)
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "encoding": "bpforms.ProteinForm",
              "structure": "A"
            },
            {
              "subunit_name": "abc_b",
              "encoding": "bpforms.ProteinForm",
              "structure": "M"
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "structure": "C[C@H]([NH3+])C(=O)O.CSCC[C@H]([NH3+])C(=O)O",
            "formula": "C8H20N2O4S",
            "mol_wt": 240.318,
            "charge": 2
        })

        # test when all structure is known (dna)
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "encoding": "bpforms.DnaForm",
              "structure": "A"
            },
            {
              "subunit_name": "abc_b",
              "encoding": "bpforms.DnaForm",
              "structure": "T"
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "structure": "OC1CC(OC1COP(=O)([O-])[O-])n1cnc2c1ncnc2N.OC1CC(OC1COP(=O)([O-])[O-])n1cc(C)c(=O)[nH]c1=O",
            "formula": "C20H25N7O14P2",
            "mol_wt": 649.402523996,
            "charge": -4
        })

        # test when all structure is known (rna)
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "encoding": "bpforms.RnaForm",
              "structure": "A"
            },
            {
              "subunit_name": "abc_b",
              "encoding": "bpforms.RnaForm",
              "structure": "U"
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "structure": "OC1C(O)C(OC1n1cnc2c1ncnc2N)COP(=O)([O-])[O-].OC1C(O)C(OC1n1ccc(=O)[nH]c1=O)COP(=O)([O-])[O-]",
            "formula": "C19H23N7O16P2",
            "mol_wt": 667.3735239959999,
            "charge": -4
        })

        # test when combination of bpform and smiles
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "encoding": "bpforms.ProteinForm",
              "structure": "A"
            },
            {
              "subunit_name": "abc_b",
              "encoding": "smiles",
              "structure": "[Zn+2]"
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "structure": "C[C@H]([NH3+])C(=O)O.[Zn+2]",
            "formula": "C3H8NO2Zn",
            "mol_wt": 155.482,
            "charge": 3
        })

        # when no structure is known, and some properties are known
        # all formula known -> formula + mol_wt
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "formula": "CH4"
            },
            {
              "subunit_name": "abc_b",
              "formula": "H2O"
            }
          ]
          })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "formula": "CH6O",
            "mol_wt": 34.058
        })

        # all mol_wt known -> mol_wt
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "mol_wt": 16
            },
            {
              "subunit_name": "abc_b",
              "mol_wt": 18
            }
          ]
          })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "mol_wt": 34.0
        })

        # all charge known -> charge
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "charge": 1
            },
            {
              "subunit_name": "abc_b",
              "charge": -1
            }
          ]
          })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "charge": 0
        })

        # when mix-and-match known information

        # some known structure + some known formula -> formula, mol_wt
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "encoding": "bpforms.ProteinForm",
              "structure": "A"
            },
            {
              "subunit_name": "abc_b",
              "formula": "CH4"
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "formula": "C4H12NO2",
            "mol_wt": 106.14500000000001
        })

        # some known structure + some known mol_wt -> mol_wt
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "encoding": "bpforms.ProteinForm",
              "structure": "A"
            },
            {
              "subunit_name": "abc_b",
              "mol_wt": 16
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "mol_wt": 106.102
        })

        # some known formula + some known mol_wt -> mol_wt
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "formula": "CH4"
            },
            {
              "subunit_name": "abc_b",
              "mol_wt": 16
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b",
            "mol_wt": 32.043
        })

        # some known formula + some known charge -> nothing
        rv = client.post('/api/bcform/', json={
          "form": "abc_a + abc_b",
          "subunits": [
            {
              "subunit_name": "abc_a",
              "formula": "CH4"
            },
            {
              "subunit_name": "abc_b",
              "charge": 1
            }
          ]
        })
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            "form": "1 * abc_a + 1 * abc_b"
        })




    # def test_get_bcform_properties_errors(self):
    #
    #     client = rest.app.test_client()
    #
    #     # invalid forms
    #     rv = client.post('/api/bcform/', json=dict(form='HELLO'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b | crosslink: [left-bond-atom: abc_c(1)-2O1 | left-displaced-atom: abc_d(1)-2H1 | right-bond-atom: abc_b(1)-3C1 | right-displaced-atom: abc_b(1)-3H1 | right-displaced-atom: abc_b(1)-3O1]'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     # invalid subunit_formulas
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_formulas='{abc_a:C5H10O, abc_b}'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_formulas='{abc_a:C5H10O}'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     # invalid subunit_mol_wts
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_mol_wts='{abc_a:86, abc_b}'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_mol_wts='{abc_a:86}'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     # invalid subunit_charges
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_charges='{abc_a:+1, abc_b:0.5}'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_charges='{abc_a:+1, abc_b}'))
    #     self.assertEqual(rv.status_code, 400)
    #
    #     rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
    #                                                 subunit_charges='{abc_a:+1}'))
    #     self.assertEqual(rv.status_code, 400)
    #
