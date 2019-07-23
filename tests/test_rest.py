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

'''
    def test_get_bcform_properties(self):
        client = rest.app.test_client()

        # test validate
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b'
        })

        # test get_formula
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_formulas='{abc_a:C5H10O, abc_b:C3H5O}'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b',
            'formula': 'C8H15O2',
        })

        # test get_mol_wt
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_mol_wts='{abc_a:86, abc_b:57}'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b',
            'mol_wt': '143.0',
        })

        # test get_charge
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_charges='{abc_a:+1, abc_b:-1}'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b',
            'charge': '0',
        })

        # # test get_formula + get_mol_wt
        # rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
        #                                             subunit_formulas='{abc_a:C5H10O, abc_b:C3H5O}',
        #                                             subunit_mol_wts='{abc_a:86, abc_b:57}'))
        # self.assertEqual(rv.status_code, 200)
        # self.assertEqual(rv.get_json(), {
        #     'form': '1 * abc_a + 1 * abc_b',
        #     'formula': 'C8H15O2',
        #     'mol_wt': '143.0',
        # })

        # test get_formula + get_charge
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_formulas='{abc_a:C5H10O, abc_b:C3H5O}',
                                                    subunit_charges='{abc_a:+1, abc_b:-1}'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b',
            'formula': 'C8H15O2',
            'charge': '0',
        })

        # test get_mol_wt + get_charge
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_mol_wts='{abc_a:86, abc_b:57}',
                                                    subunit_charges='{abc_a:+1, abc_b:-1}'))
        self.assertEqual(rv.status_code, 200)
        self.assertEqual(rv.get_json(), {
            'form': '1 * abc_a + 1 * abc_b',
            'mol_wt': '143.0',
            'charge': '0',
        })

        # # test get_formula + get_mol_wt + get_charge
        # rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
        #                                             subunit_formulas='{abc_a:C5H10O, abc_b:C3H5O}',
        #                                             subunit_mol_wts='{abc_a:86, abc_b:57}',
        #                                             subunit_charges='{abc_a:+1, abc_b:-1}'))
        # self.assertEqual(rv.status_code, 200)
        # self.assertEqual(rv.get_json(), {
        #     'form': '1 * abc_a + 1 * abc_b',
        #     'formula': 'C8H15O2',
        #     'mol_wt': '143.0',
        #     'charge': '0',
        # })

    def test_get_bcform_properties_errors(self):

        client = rest.app.test_client()

        # invalid forms
        rv = client.post('/api/bcform/', json=dict(form='HELLO'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b | crosslink: [left-bond-atom: abc_c(1)-2O1 | left-displaced-atom: abc_d(1)-2H1 | right-bond-atom: abc_b(1)-3C1 | right-displaced-atom: abc_b(1)-3H1 | right-displaced-atom: abc_b(1)-3O1]'))
        self.assertEqual(rv.status_code, 400)

        # invalid subunit_formulas
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_formulas='{abc_a:C5H10O, abc_b}'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_formulas='{abc_a:C5H10O}'))
        self.assertEqual(rv.status_code, 400)

        # invalid subunit_mol_wts
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_mol_wts='{abc_a:86, abc_b}'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_mol_wts='{abc_a:86}'))
        self.assertEqual(rv.status_code, 400)

        # invalid subunit_charges
        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_charges='{abc_a:+1, abc_b:0.5}'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_charges='{abc_a:+1, abc_b}'))
        self.assertEqual(rv.status_code, 400)

        rv = client.post('/api/bcform/', json=dict(form='abc_a + abc_b',
                                                    subunit_charges='{abc_a:+1}'))
        self.assertEqual(rv.status_code, 400)
'''
