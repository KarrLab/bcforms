""" Test of bcforms.core

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-6-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bcforms import core
import unittest


class BcFormTestCase(unittest.TestCase):

    def test_init(self):
        bc_form = core.BcForm()
        self.assertEqual(bc_form.subunits, [])
        self.assertEqual(bc_form.crosslinks, [])

    def test_set_subunits(self):
        bc_form = core.BcForm()

        bc_form.subunits = []

        with self.assertRaises(ValueError):
            bc_form.subunits = None


    def test_set_crosslinks(self):

        bc_form = core.BcForm()

        bc_form.crosslinks = []

        with self.assertRaises(ValueError):
            bc_form.crosslinks = None
