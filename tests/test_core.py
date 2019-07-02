""" Test of bcforms.core

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-6-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bcforms import core
import unittest

class AtomTestCase(unittest.TestCase):

    def test_init(self):
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        self.assertEqual(atom.subunit, 'abc')
        self.assertEqual(atom.subunit_idx, 1)
        self.assertEqual(atom.element, 'H')
        self.assertEqual(atom.position, 1)
        self.assertEqual(atom.monomer, 10)
        self.assertEqual(atom.charge, 0)

    def test_set_subunit(self):
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.subunit = 'def'
        with self.assertRaises(ValueError):
            atom.subunit = None

    def test_set_subunit_idx(self):
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.subunit_idx = 2
        with self.assertRaises(ValueError):
            atom.subunit_idx = -1
        with self.assertRaises(ValueError):
            atom.subunit_idx = None

    def test_set_element(self):
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.element = 'N'
        with self.assertRaises(ValueError):
            atom.element = None

    def test_set_position(self):
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.position = 2
        with self.assertRaises(ValueError):
            atom.position = -1
        with self.assertRaises(ValueError):
            atom.position = None

    def test_set_monomer(self):
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.monomer = 2
        with self.assertRaises(ValueError):
            atom.monomer = -1
        with self.assertRaises(ValueError):
            atom.monomer = None

    def test_set_charge(self):

        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.charge = 2
        with self.assertRaises(ValueError):
            atom.charge = None

class CrosslinkTestCase(unittest.TestCase):

    def test_init(self):
        crosslink = core.Crosslink()
        self.assertEqual(crosslink.left_bond_atoms, [])
        self.assertEqual(crosslink.right_bond_atoms, [])
        self.assertEqual(crosslink.left_displaced_atoms, [])
        self.assertEqual(crosslink.right_displaced_atoms, [])

    def test_set_left_bond_atoms(self):
        crosslink = core.Crosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.left_bond_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.left_bond_atoms = None

    def test_set_right_bond_atoms(self):
        crosslink = core.Crosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.right_bond_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.right_bond_atoms = None

    def test_set_left_displaced_atoms(self):
        crosslink = core.Crosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.left_displaced_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.left_displaced_atoms = None

    def test_set_right_displaced_atoms(self):
        crosslink = core.Crosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.right_displaced_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.right_displaced_atoms = None



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

    def test_from_str(self):

        bc_form_1 = core.BcForm().from_str('2 * abc_a + 3 * abc_b')

        self.assertEqual(len(bc_form_1.subunits), 2)
        self.assertEqual(bc_form_1.subunits[0]['id'], 'abc_a')
        self.assertEqual(bc_form_1.subunits[0]['stoichiometry'], 2)
        self.assertEqual(bc_form_1.subunits[1]['id'], 'abc_b')
        self.assertEqual(bc_form_1.subunits[1]['stoichiometry'], 3)
        self.assertEqual(bc_form_1.crosslinks, [])

        bc_form_2 = core.BcForm().from_str('bmp2_a + bmp2_a | crosslink: [left-bond-atom: bmp2_a(1)-362S1 | left-displaced-atom: bmp2_a(1)-362H1 | right-bond-atom: bmp2_a(2)-362S1 | right-displaced-atom: bmp2_a(2)-362H1]')

        self.assertEqual(len(bc_form_2.crosslinks), 1)
        self.assertEqual(len(bc_form_2.crosslinks[0].left_bond_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].left_bond_atoms[0].element, 'S')
        self.assertEqual(len(bc_form_2.crosslinks[0].right_bond_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].right_bond_atoms[0].element, 'S')
        self.assertEqual(len(bc_form_2.crosslinks[0].left_displaced_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].left_displaced_atoms[0].element, 'H')
        self.assertEqual(len(bc_form_2.crosslinks[0].right_displaced_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].right_displaced_atoms[0].element, 'H')

    def test_from_set(self):

        bc_form = core.BcForm().from_set([{'id': 'abc_a', 'stoichiometry': 2}, {'id': 'abc_b', 'stoichiometry': 3}])

        self.assertEqual(len(bc_form.subunits), 2)
        self.assertEqual(bc_form.subunits[0]['id'], 'abc_a')
        self.assertEqual(bc_form.subunits[0]['stoichiometry'], 2)
        self.assertEqual(bc_form.subunits[1]['id'], 'abc_b')
        self.assertEqual(bc_form.subunits[1]['stoichiometry'], 3)
        self.assertEqual(bc_form.crosslinks, [])

    def test_clean(self):

        bc_form = core.BcForm().from_str('abc + abc + abc')
        self.assertEqual(len(bc_form.subunits), 1)
        self.assertEqual(bc_form.subunits[0]['id'], 'abc')
        self.assertEqual(bc_form.subunits[0]['stoichiometry'], 3)
