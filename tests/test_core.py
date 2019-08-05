""" Test of bcforms.core

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-6-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bcforms import core
import openbabel
import unittest
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
import bpforms
import bpforms.core
import bpforms.alphabet.protein

class SubunitTestCase(unittest.TestCase):

    def test_init(self):
        subunit_1 = core.Subunit(id='abc', stoichiometry=2)
        self.assertEqual(subunit_1.id, 'abc')
        self.assertEqual(subunit_1.stoichiometry, 2)
        self.assertIsNone(subunit_1.structure)
        self.assertIsNone(subunit_1.formula)
        self.assertIsNone(subunit_1.mol_wt)
        self.assertIsNone(subunit_1.charge)

        subunit_2 = core.Subunit(id='abc', stoichiometry=2, structure=bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(subunit_2.id, 'abc')
        self.assertEqual(subunit_2.stoichiometry, 2)
        self.assertEqual(len(subunit_2.structure), 2)

    def test_set_id(self):
        subunit = core.Subunit(id='abc', stoichiometry=2)
        subunit.id = 'def'
        self.assertEqual(subunit.id, 'def')
        with self.assertRaises(ValueError):
            subunit.id = None

    def test_set_stoichiometry(self):
        subunit = core.Subunit(id='abc', stoichiometry=3)
        subunit.stoichiometry = 3
        self.assertEqual(subunit.stoichiometry, 3)
        with self.assertRaises(ValueError):
            subunit.stoichiometry = None

    def test_set_structure(self):
        subunit = core.Subunit(id='abc', stoichiometry=3)
        subunit.structure = bpforms.ProteinForm().from_str('AA')
        self.assertIsNotNone(subunit.formula)
        self.assertIsNotNone(subunit.mol_wt)
        self.assertIsNotNone(subunit.charge)
        subunit.structure = openbabel.OBMol()
        subunit.structure = None
        with self.assertRaises(ValueError):
            subunit.structure = 123
        with self.assertRaises(ValueError):
            subunit.structure = 'CH3'

    def test_set_formula(self):
        subunit = core.Subunit(id='abc', stoichiometry=3)
        subunit.formula = EmpiricalFormula('CH4')
        self.assertIsNotNone(subunit.mol_wt)
        subunit.formula = None
        with self.assertRaises(ValueError):
            subunit.formula = 123
        with self.assertRaises(ValueError):
            subunit.structure = bpforms.ProteinForm().from_str('AA')
            subunit.formula = EmpiricalFormula('CH4')

    def test_set_mol_wt(self):
        subunit = core.Subunit(id='abc', stoichiometry=3)
        subunit.mol_wt = 16.0
        subunit.mol_wt = 12
        subunit.mol_wt = None
        with self.assertRaises(ValueError):
            subunit.mol_wt = 'string'
        with self.assertRaises(ValueError):
            subunit.formula = EmpiricalFormula('CH4')
            subunit.mol_wt = 12.0
        with self.assertRaises(ValueError):
            subunit.mol_wt = -5

    def test_set_charge(self):
        subunit = core.Subunit(id='abc', stoichiometry=3)
        subunit.charge = 2
        subunit.charge = None
        with self.assertRaises(ValueError):
            subunit.charge = 0.5
        with self.assertRaises(ValueError):
            subunit.structure = bpforms.ProteinForm().from_str('AA')
            subunit.charge = -3

    def test_str(self):
        subunit = core.Subunit(id='abc', stoichiometry=3)
        self.assertEqual(str(subunit), '3 * abc')

    def test_equal(self):

        subunit_1 = core.Subunit(id='abc', stoichiometry=3)
        subunit_2 = core.Subunit(id='abc', stoichiometry=3)
        subunit_3 = core.Subunit(id='def', stoichiometry=3)
        subunit_4 = core.Subunit(id='abc', stoichiometry=2)

        self.assertFalse(subunit_1.is_equal('abc'))
        self.assertTrue(subunit_1.is_equal(subunit_1))
        self.assertTrue(subunit_1.is_equal(subunit_2))
        self.assertFalse(subunit_1.is_equal(subunit_3))
        self.assertFalse(subunit_1.is_equal(subunit_4))

    def test_get_formula(self):
        subunit_1 = core.Subunit(id='abc', stoichiometry=2)
        self.assertEqual(subunit_1.get_formula(formula=EmpiricalFormula('CH4')), EmpiricalFormula('C2H8'))

        subunit_2 = core.Subunit(id='abc', stoichiometry=2)
        self.assertIsNone(subunit_2.get_formula())

        mol = openbabel.OBMol()
        a = mol.NewAtom()
        a.SetAtomicNum(12)
        subunit_3 = core.Subunit(id='mg', stoichiometry=1, structure=mol)
        self.assertEqual(subunit_3.get_formula(), EmpiricalFormula('Mg'))

        subunit_4 = core.Subunit(id='aa', stoichiometry=1, structure=bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(subunit_4.get_formula(), EmpiricalFormula('C6H13N2O3'))

    def test_get_mol_wt(self):
        subunit_1 = core.Subunit(id='abc', stoichiometry=2)
        self.assertEqual(subunit_1.get_mol_wt(mol_wt=32.0), 64.0)

        subunit_2 = core.Subunit(id='abc', stoichiometry=2)
        self.assertIsNone(subunit_2.get_mol_wt())

        mol = openbabel.OBMol()
        a = mol.NewAtom()
        a.SetAtomicNum(12)
        subunit_3 = core.Subunit(id='mg', stoichiometry=1, structure=mol)
        self.assertAlmostEqual(subunit_3.get_mol_wt(), 24, places=0)

        subunit_4 = core.Subunit(id='aa', stoichiometry=2, structure=bpforms.ProteinForm().from_str('AA'))
        self.assertAlmostEqual(subunit_4.get_mol_wt(), 322.362, places=3)

    def test_get_charge(self):
        subunit_1 = core.Subunit(id='abc', stoichiometry=2)
        self.assertEqual(subunit_1.get_charge(charge=1), 2)

        subunit_2 = core.Subunit(id='abc', stoichiometry=2)
        self.assertIsNone(subunit_2.get_charge())

        mol = openbabel.OBMol()
        a = mol.NewAtom()
        a.SetAtomicNum(12)
        subunit_3 = core.Subunit(id='mg', stoichiometry=2, structure=mol)
        self.assertEqual(subunit_3.get_charge(), 0)

        subunit_4 = core.Subunit(id='aa', stoichiometry=1, structure=bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(subunit_4.get_charge(), 1)

    def test_get_structure(self):

        subunit_1 = core.Subunit(id='aa', stoichiometry=2, structure=bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(OpenBabelUtils.export(subunit_1.get_structure()[0], 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)O.C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)O')
        self.assertTrue(isinstance(subunit_1.get_structure()[1][2], dict))

        ob_mol = openbabel.OBMol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('smi')
        conversion.ReadString(ob_mol, 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)[O-]')
        subunit_2 = core.Subunit(id='aa', stoichiometry=1, structure=ob_mol)
        self.assertEqual(OpenBabelUtils.export(subunit_2.get_structure()[0], 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)[O-]')
        self.assertEqual(len(subunit_2.get_structure()[1][1][1]['monomer']),16)

        subunit_3 = core.Subunit(id='aa', stoichiometry=1)
        with self.assertRaises(ValueError):
            subunit_3.get_structure()

        subunit_4 = core.Subunit(id='aa', stoichiometry=1, structure='C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)[O-]')
        self.assertEqual(OpenBabelUtils.export(subunit_4.get_structure()[0], 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)[O-]')

    def test_export(self):
        subunit_1 = core.Subunit(id='aa', stoichiometry=2, structure=bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(subunit_1.export(), 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)O.C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)O')

class AtomTestCase(unittest.TestCase):

    def test_init(self):
        atom_1 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0, component_type='monomer')
        self.assertEqual(atom_1.subunit, 'abc')
        self.assertEqual(atom_1.subunit_idx, 1)
        self.assertEqual(atom_1.element, 'H')
        self.assertEqual(atom_1.position, 1)
        self.assertEqual(atom_1.monomer, 10)
        self.assertEqual(atom_1.charge, 0)
        self.assertEqual(atom_1.component_type, 'monomer')

        atom_2 = core.Atom(subunit='abc', subunit_idx=None, element='H', position=1, monomer=10, charge=0)
        self.assertIsNone(atom_2.subunit_idx)


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

    def test_set_component_type(self):

        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom.component_type = 'backbone'
        with self.assertRaises(ValueError):
            atom.component_type = 'm'

    def test_str(self):
        atom_1 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        self.assertEqual(str(atom_1), 'abc(1)-10H1')
        atom_2 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=1)
        self.assertEqual(str(atom_2), 'abc(1)-10H1+1')

    def test_is_equal(self):

        atom_1 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom_2 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom_3 = core.Atom(subunit='def', subunit_idx=1, element='H', position=1, monomer=10, charge=0)

        self.assertTrue(atom_1.is_equal(atom_1))
        self.assertFalse(atom_1.is_equal('atom'))
        self.assertTrue(atom_1.is_equal(atom_2))
        self.assertFalse(atom_1.is_equal(atom_3))

        atom_4 = core.Atom(subunit='abc', subunit_idx=3, element='H', position=1, monomer=10, charge=0)
        atom_5 = core.Atom(subunit='abc', subunit_idx=None, element='H', position=1, monomer=10, charge=0)
        self.assertFalse(atom_4.is_equal(atom_5))

class InlineCrosslinkTestCase(unittest.TestCase):

    def test_init(self):
        crosslink = core.InlineCrosslink()
        self.assertEqual(crosslink.l_bond_atoms, [])
        self.assertEqual(crosslink.r_bond_atoms, [])
        self.assertEqual(crosslink.l_displaced_atoms, [])
        self.assertEqual(crosslink.r_displaced_atoms, [])

        crosslink = core.InlineCrosslink([], [], [], [])
        self.assertEqual(crosslink.l_bond_atoms, [])
        self.assertEqual(crosslink.r_bond_atoms, [])
        self.assertEqual(crosslink.l_displaced_atoms, [])
        self.assertEqual(crosslink.r_displaced_atoms, [])

    def test_set_l_bond_atoms(self):
        crosslink = core.InlineCrosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.l_bond_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.l_bond_atoms = None

    def test_set_r_bond_atoms(self):
        crosslink = core.InlineCrosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.r_bond_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.r_bond_atoms = None

    def test_set_l_displaced_atoms(self):
        crosslink = core.InlineCrosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.l_displaced_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.l_displaced_atoms = None

    def test_set_r_displaced_atoms(self):
        crosslink = core.InlineCrosslink()
        atom = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.r_displaced_atoms.append(atom)
        with self.assertRaises(ValueError):
            crosslink.r_displaced_atoms = None

    def test_str(self):
        crosslink = core.InlineCrosslink()
        atom_1 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.l_bond_atoms.append(atom_1)
        atom_2 = core.Atom(subunit='def', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        crosslink.r_bond_atoms.append(atom_2)

        self.assertEqual(str(crosslink), 'x-link: [ l-bond-atom: abc(1)-10H1 | r-bond-atom: def(1)-10H1 ]')

    def test_is_equal(self):
        atom_1 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom_2 = core.Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
        atom_3 = core.Atom(subunit='def', subunit_idx=1, element='H', position=1, monomer=10, charge=0)

        crosslink_1 = core.InlineCrosslink()
        crosslink_1.l_bond_atoms.append(atom_1)
        crosslink_1.r_bond_atoms.append(atom_2)

        crosslink_2 = core.InlineCrosslink()
        crosslink_2.l_bond_atoms.append(atom_1)
        crosslink_2.r_bond_atoms.append(atom_2)

        crosslink_3 = core.InlineCrosslink()
        crosslink_3.l_bond_atoms.append(atom_1)
        crosslink_3.r_bond_atoms.append(atom_2)
        crosslink_3.r_bond_atoms.append(atom_3)

        crosslink_4 = core.InlineCrosslink()
        crosslink_4.l_bond_atoms.append(atom_1)
        crosslink_4.r_bond_atoms.append(atom_3)

        self.assertTrue(crosslink_1.is_equal(crosslink_1))
        self.assertFalse(crosslink_1.is_equal(atom_1))
        self.assertTrue(crosslink_1.is_equal(crosslink_2))
        self.assertFalse(crosslink_1.is_equal(crosslink_3))
        self.assertFalse(crosslink_1.is_equal(crosslink_4))

        xlink_5_a = core.InlineCrosslink()
        xlink_5_a.l_bond_atoms.append(core.Atom(subunit='sub_1', subunit_idx=None, element='S', position=11, monomer=1, charge=0))
        xlink_5_a.l_displaced_atoms.append(core.Atom(subunit='sub_1', subunit_idx=None, element='H', position=11, monomer=1, charge=0))
        xlink_5_a.r_bond_atoms.append(core.Atom(subunit='sub_2', subunit_idx=None, element='S', position=11, monomer=1, charge=0))
        xlink_5_a.r_displaced_atoms.append(core.Atom(subunit='sub_2', subunit_idx=None, element='H', position=11, monomer=1, charge=0))

        xlink_5_b = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=1, r_monomer=1)
        self.assertTrue(xlink_5_a.is_equal(xlink_5_b))

class AbstractedCrosslinkTestCase(unittest.TestCase):

    def test_init(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        self.assertEqual(xlink_1.type, 'disulfide')
        self.assertEqual(xlink_1.l_subunit, 'sub_1')
        self.assertIsNone(xlink_1.l_subunit_idx)
        self.assertEqual(xlink_1.l_monomer, 3)
        self.assertEqual(xlink_1.r_subunit, 'sub_2')
        self.assertIsNone(xlink_1.r_subunit_idx)
        self.assertEqual(xlink_1.r_monomer, 1)

        xlink_2 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', l_subunit_idx=2, r_subunit='sub_2', r_subunit_idx=1, l_monomer=3, r_monomer=1)
        self.assertEqual(xlink_2.l_subunit_idx, 2)
        self.assertEqual(xlink_2.r_subunit_idx, 1)

    def test_set_type(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.type = 'glycyl_lysine_isopeptide'
        self.assertEqual(xlink_1.type, 'glycyl_lysine_isopeptide')
        with self.assertRaises(ValueError):
            xlink_1.type = None

    def test_set_l_subunit(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.l_subunit = 'sub_3'
        self.assertEqual(xlink_1.l_subunit, 'sub_3')
        with self.assertRaises(ValueError):
            xlink_1.l_subunit = None

    def test_set_l_subunit_idx(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.l_subunit_idx = 2
        self.assertEqual(xlink_1.l_subunit_idx, 2)
        with self.assertRaises(ValueError):
            xlink_1.l_subunit_idx = 'AA'
        with self.assertRaises(ValueError):
            xlink_1.l_subunit_idx = -1

    def test_set_l_monomer(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.l_monomer = 2
        self.assertEqual(xlink_1.l_monomer, 2)
        with self.assertRaises(ValueError):
            xlink_1.l_monomer = None

    def test_set_r_subunit(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.r_subunit = 'sub_3'
        self.assertEqual(xlink_1.r_subunit, 'sub_3')
        with self.assertRaises(ValueError):
            xlink_1.r_subunit = None

    def test_set_r_subunit_idx(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.r_subunit_idx = 2
        self.assertEqual(xlink_1.r_subunit_idx, 2)
        with self.assertRaises(ValueError):
            xlink_1.r_subunit_idx = 'AA'
        with self.assertRaises(ValueError):
            xlink_1.r_subunit_idx = -1

    def test_set_r_monomer(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        xlink_1.r_monomer = 2
        self.assertEqual(xlink_1.r_monomer, 2)
        with self.assertRaises(ValueError):
            xlink_1.r_monomer = None

    def test_str(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=3, r_monomer=1)
        self.assertEqual(str(xlink_1), 'x-link: [ type: disulfide | l-monomer: sub_1-3 | r-monomer: sub_2-1 ]')

        xlink_2 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', l_subunit_idx=2, r_subunit='sub_2', r_subunit_idx=1, l_monomer=3, r_monomer=1)
        self.assertEqual(str(xlink_2), 'x-link: [ type: disulfide | l-monomer: sub_1(2)-3 | r-monomer: sub_2(1)-1 ]')

    def test_is_equal(self):
        xlink_1_a = core.InlineCrosslink()
        xlink_1_a.l_bond_atoms.append(core.Atom(subunit='sub_1', subunit_idx=None, element='S', position=11, monomer=1, charge=0))
        xlink_1_a.l_displaced_atoms.append(core.Atom(subunit='sub_1', subunit_idx=None, element='H', position=11, monomer=1, charge=0))
        xlink_1_a.r_bond_atoms.append(core.Atom(subunit='sub_2', subunit_idx=None, element='S', position=11, monomer=1, charge=0))
        xlink_1_a.r_displaced_atoms.append(core.Atom(subunit='sub_2', subunit_idx=None, element='H', position=11, monomer=1, charge=0))

        xlink_1_b = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=1, r_monomer=1)
        self.assertTrue(xlink_1_b.is_equal(xlink_1_a))

        xlink_1_c = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', r_subunit='sub_2', l_monomer=1, r_monomer=1)
        self.assertTrue(xlink_1_b.is_equal(xlink_1_c))

        xlink_2 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_3', r_subunit='sub_2', l_monomer=1, r_monomer=1)
        self.assertFalse(xlink_1_b.is_equal(xlink_2))

    def test_get_l_bond_atoms(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', l_subunit_idx=1, r_subunit='sub_1', r_subunit_idx=2, l_monomer=1, r_monomer=1)
        self.assertEqual(len(xlink_1.get_l_bond_atoms()), 1)
        self.assertEqual(xlink_1.get_l_bond_atoms()[0].subunit, 'sub_1')
        self.assertEqual(xlink_1.get_l_bond_atoms()[0].subunit_idx, 1)
        self.assertEqual(xlink_1.get_l_bond_atoms()[0].monomer, 1)
        self.assertEqual(xlink_1.get_l_bond_atoms()[0].element, 'S')
        self.assertEqual(xlink_1.get_l_bond_atoms()[0].position, 11)
        self.assertEqual(xlink_1.get_l_bond_atoms()[0].charge, 0)

    def test_get_r_bond_atoms(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', l_subunit_idx=1, r_subunit='sub_1', r_subunit_idx=2, l_monomer=1, r_monomer=1)
        self.assertEqual(len(xlink_1.get_r_bond_atoms()), 1)
        self.assertEqual(xlink_1.get_r_bond_atoms()[0].subunit, 'sub_1')
        self.assertEqual(xlink_1.get_r_bond_atoms()[0].subunit_idx, 2)
        self.assertEqual(xlink_1.get_r_bond_atoms()[0].monomer, 1)
        self.assertEqual(xlink_1.get_r_bond_atoms()[0].element, 'S')
        self.assertEqual(xlink_1.get_r_bond_atoms()[0].position, 11)
        self.assertEqual(xlink_1.get_r_bond_atoms()[0].charge, 0)

    def test_get_l_displaced_atoms(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', l_subunit_idx=1, r_subunit='sub_1', r_subunit_idx=2, l_monomer=1, r_monomer=1)
        self.assertEqual(len(xlink_1.get_l_displaced_atoms()), 1)
        self.assertEqual(xlink_1.get_l_displaced_atoms()[0].subunit, 'sub_1')
        self.assertEqual(xlink_1.get_l_displaced_atoms()[0].subunit_idx, 1)
        self.assertEqual(xlink_1.get_l_displaced_atoms()[0].monomer, 1)
        self.assertEqual(xlink_1.get_l_displaced_atoms()[0].element, 'H')
        self.assertEqual(xlink_1.get_l_displaced_atoms()[0].position, 11)
        self.assertEqual(xlink_1.get_l_displaced_atoms()[0].charge, 0)

    def test_get_r_displaced_atoms(self):
        xlink_1 = core.AbstractedCrosslink(type='disulfide', l_subunit='sub_1', l_subunit_idx=1, r_subunit='sub_1', r_subunit_idx=2, l_monomer=1, r_monomer=1)
        self.assertEqual(len(xlink_1.get_r_displaced_atoms()), 1)
        self.assertEqual(xlink_1.get_r_displaced_atoms()[0].subunit, 'sub_1')
        self.assertEqual(xlink_1.get_r_displaced_atoms()[0].subunit_idx, 2)
        self.assertEqual(xlink_1.get_r_displaced_atoms()[0].monomer, 1)
        self.assertEqual(xlink_1.get_r_displaced_atoms()[0].element, 'H')
        self.assertEqual(xlink_1.get_r_displaced_atoms()[0].position, 11)
        self.assertEqual(xlink_1.get_r_displaced_atoms()[0].charge, 0)

class BcFormTestCase(unittest.TestCase):

    def test_init(self):
        bc_form = core.BcForm()
        self.assertEqual(bc_form.subunits, [])
        self.assertEqual(bc_form.crosslinks, [])

        bc_form = core.BcForm([],[])
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

    def test_str(self):

        s1 = '2 * abc_a + 3 * abc_b'
        bc_form_1 = core.BcForm().from_str(s1)
        self.assertEqual(s1, str(bc_form_1))

        s2 = '2 * bmp2_a | x-link: [ l-bond-atom: bmp2_a(1)-362S1 | l-displaced-atom: bmp2_a(1)-362H1 | r-bond-atom: bmp2_a(2)-362S1 | r-displaced-atom: bmp2_a(2)-362H1 ]'
        bc_form_2 = core.BcForm().from_str(s2)
        self.assertEqual(s2, str(bc_form_2))

        bc_form_3 = core.BcForm()
        bc_form_3.subunits.append(core.Subunit(id='bmp2_a', stoichiometry=1))
        bc_form_3.crosslinks.append(core.InlineCrosslink(l_bond_atoms=[
            core.Atom(subunit='bmp2_a', subunit_idx=None, element='H', position=1, monomer=10, charge=0)]))
        self.assertEqual(str(bc_form_3), '1 * bmp2_a | x-link: [ l-bond-atom: bmp2_a-10H1 ]')

        bc_form_4 = core.BcForm().from_str('unit_1 + unit_2'
                                          '| x-link: [ type: glycyl_lysine_isopeptide |'
                                                     ' l-monomer: unit_1-1 |'
                                                     ' r-monomer: unit_2-2 ]')
        self.assertEqual(str(bc_form_4), '1 * unit_1 + 1 * unit_2 | x-link: [ type: glycyl_lysine_isopeptide | l-monomer: unit_1-1 | r-monomer: unit_2-2 ]')

    def test_from_str(self):

        bc_form_1 = core.BcForm().from_str('2 * abc_a + 3 * abc_b')

        self.assertEqual(len(bc_form_1.subunits), 2)
        self.assertEqual(bc_form_1.subunits[0].id, 'abc_a')
        self.assertEqual(bc_form_1.subunits[0].stoichiometry, 2)
        self.assertEqual(bc_form_1.subunits[1].id, 'abc_b')
        self.assertEqual(bc_form_1.subunits[1].stoichiometry, 3)
        self.assertEqual(bc_form_1.crosslinks, [])

        bc_form_2 = core.BcForm().from_str('bmp2_a + bmp2_a | x-link: [l-bond-atom: bmp2_a(1)-362S1 | l-displaced-atom: bmp2_a(1)-362H1 | r-bond-atom: bmp2_a(2)-362S1 | r-displaced-atom: bmp2_a(2)-362H1]')

        self.assertEqual(len(bc_form_2.crosslinks), 1)
        self.assertEqual(len(bc_form_2.crosslinks[0].l_bond_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].l_bond_atoms[0].element, 'S')
        self.assertEqual(len(bc_form_2.crosslinks[0].r_bond_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].r_bond_atoms[0].element, 'S')
        self.assertEqual(len(bc_form_2.crosslinks[0].l_displaced_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].l_displaced_atoms[0].element, 'H')
        self.assertEqual(len(bc_form_2.crosslinks[0].r_displaced_atoms), 1)
        self.assertEqual(bc_form_2.crosslinks[0].r_displaced_atoms[0].element, 'H')

        bc_form_3 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a-2O1 | l-displaced-atom: abc_a-2H1 | r-bond-atom: abc_b-3C1 | r-displaced-atom: abc_b-3H1 | r-displaced-atom: abc_b-3O1]')
        self.assertIsNone(bc_form_3.crosslinks[0].l_bond_atoms[0].subunit_idx)
        self.assertEqual(bc_form_3.crosslinks[0].l_bond_atoms[0].charge, 0)

        bc_form_4 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a-2O1+1 | l-displaced-atom: abc_a-2H1 | r-bond-atom: abc_b-3C1 | r-displaced-atom: abc_b-3H1 | r-displaced-atom: abc_b-3O1]')
        self.assertIsNone(bc_form_4.crosslinks[0].l_bond_atoms[0].subunit_idx)
        self.assertEqual(len(bc_form_4.crosslinks[0].l_bond_atoms), 1)
        self.assertEqual(bc_form_4.crosslinks[0].l_bond_atoms[0].element, 'O')
        self.assertEqual(bc_form_4.crosslinks[0].l_bond_atoms[0].charge, 1)

        bc_form_5 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a-2O1m+1 | l-displaced-atom: abc_a-2H1m | r-bond-atom: abc_b-3C1 | r-displaced-atom: abc_b-3H1+1 | r-displaced-atom: abc_b-3O1b-1]')
        self.assertEqual(bc_form_5.crosslinks[0].l_bond_atoms[0].component_type, 'monomer')
        self.assertEqual(bc_form_5.crosslinks[0].l_bond_atoms[0].charge, 1)
        self.assertEqual(bc_form_5.crosslinks[0].l_displaced_atoms[0].component_type, 'monomer')
        self.assertEqual(bc_form_5.crosslinks[0].l_displaced_atoms[0].charge, 0)
        self.assertEqual(bc_form_5.crosslinks[0].r_bond_atoms[0].component_type, 'monomer')
        self.assertEqual(bc_form_5.crosslinks[0].r_bond_atoms[0].charge, 0)
        self.assertEqual(bc_form_5.crosslinks[0].r_displaced_atoms[0].component_type, 'monomer')
        self.assertEqual(bc_form_5.crosslinks[0].r_displaced_atoms[0].charge, 1)
        self.assertEqual(bc_form_5.crosslinks[0].r_displaced_atoms[1].component_type, 'backbone')
        self.assertEqual(bc_form_5.crosslinks[0].r_displaced_atoms[1].charge, -1)

        bc_form_6_a = core.BcForm().from_str('2 * unit_2'
                                          '| x-link: [ l-bond-atom: unit_2(1)-1S11 |'
                                                     ' l-displaced-atom: unit_2(1)-1H11 |'
                                                     ' r-bond-atom: unit_2(2)-1S11 |'
                                                     ' r-displaced-atom: unit_2(2)-1H11 ]')
        bc_form_6_a.set_subunit_attribute('unit_2', 'structure', bpforms.ProteinForm().from_str('C'))
        bc_form_6_b = core.BcForm().from_str('2 * unit_2'
                                          '| x-link: [ type: disulfide |'
                                                     ' l-monomer: unit_2(1)-1 |'
                                                     ' r-monomer: unit_2(2)-1 ]')
        bc_form_6_b.set_subunit_attribute('unit_2', 'structure', bpforms.ProteinForm().from_str('C'))
        self.assertEqual(bc_form_6_a.export(), bc_form_6_b.export())

        bc_form_7_a = core.BcForm().from_str('2 * unit_2'
                                          '| x-link: [ l-bond-atom: unit_2(1)-2S11 |'
                                                     ' l-displaced-atom: unit_2(1)-2H11 |'
                                                     ' r-bond-atom: unit_2(2)-4S11 |'
                                                     ' r-displaced-atom: unit_2(2)-4H11 ]'
                                          '| x-link: [ l-bond-atom: unit_2(1)-4S11 |'
                                                     ' l-displaced-atom: unit_2(1)-4H11 |'
                                                     ' r-bond-atom: unit_2(2)-2S11 |'
                                                     ' r-displaced-atom: unit_2(2)-2H11 ]')
        bc_form_7_a.set_subunit_attribute('unit_2', 'structure', bpforms.ProteinForm().from_str('ACAC'))
        bc_form_7_b = core.BcForm().from_str('2 * unit_2'
                                          '| x-link: [ type: disulfide |'
                                                     ' l-monomer: unit_2(1)-2 |'
                                                     ' r-monomer: unit_2(2)-4 ]'
                                          '| x-link: [ type: disulfide |'
                                                     ' l-monomer: unit_2(1)-4 |'
                                                     ' r-monomer: unit_2(2)-2 ]')
        bc_form_7_b.set_subunit_attribute('unit_2', 'structure', bpforms.ProteinForm().from_str('ACAC'))
        self.assertEqual(bc_form_7_a.export(), bc_form_7_b.export())

        bc_form_8_a = core.BcForm().from_str('unit_1 + unit_2'
                                          '| x-link: [ l-bond-atom: unit_1-1C2 |'
                                                     ' r-bond-atom: unit_2-2N1-1 |'
                                                     ' l-displaced-atom: unit_1-1O1 |'
                                                     ' l-displaced-atom: unit_1-1H1 |'
                                                     ' r-displaced-atom: unit_2-2H1+1 |'
                                                     ' r-displaced-atom: unit_2-2H1 ]')
        bc_form_8_a.set_subunit_attribute('unit_1', 'structure', bpforms.ProteinForm().from_str('G'))
        bc_form_8_a.set_subunit_attribute('unit_2', 'structure', bpforms.ProteinForm().from_str('CKA'))
        bc_form_8_b = core.BcForm().from_str('unit_1 + unit_2'
                                          '| x-link: [ type: glycyl_lysine_isopeptide |'
                                                     ' l-monomer: unit_1-1 |'
                                                     ' r-monomer: unit_2-2 ]')
        bc_form_8_b.set_subunit_attribute('unit_1', 'structure', bpforms.ProteinForm().from_str('G'))
        bc_form_8_b.set_subunit_attribute('unit_2', 'structure', bpforms.ProteinForm().from_str('CKA'))
        self.assertEqual(bc_form_8_a.export(), bc_form_8_b.export())


    def test_from_set(self):

        bc_form = core.BcForm().from_set([{'id': 'abc_a', 'stoichiometry': 2}, {'id': 'abc_b', 'stoichiometry': 3}])

        self.assertEqual(len(bc_form.subunits), 2)
        self.assertEqual(bc_form.subunits[0].id, 'abc_a')
        self.assertEqual(bc_form.subunits[0].stoichiometry, 2)
        self.assertEqual(bc_form.subunits[1].id, 'abc_b')
        self.assertEqual(bc_form.subunits[1].stoichiometry, 3)
        self.assertEqual(bc_form.crosslinks, [])

        with self.assertRaises(ValueError):
            bad_form = core.BcForm().from_set([{'stoichiometry':2}, {'id': 'abc_b', 'stoichiometry': 3}])
        with self.assertRaises(ValueError):
            bad_form = core.BcForm().from_set([{'id':'abc'}, {'id': 'abc_b', 'stoichiometry': 3}])


    def test_clean(self):

        bc_form = core.BcForm().from_str('abc + abc + abc')
        bc_form_sub = core.BcForm().from_str('def + def')
        bc_form.subunits.append(bc_form_sub)
        bc_form.clean()
        self.assertEqual(len(bc_form.subunits), 2)
        self.assertEqual(bc_form.subunits[0].id, 'abc')
        self.assertEqual(bc_form.subunits[0].stoichiometry, 3)
        self.assertEqual(len(bc_form.subunits[1].subunits), 1)
        self.assertEqual(bc_form.subunits[1].subunits[0].id, 'def')
        self.assertEqual(bc_form.subunits[1].subunits[0].stoichiometry, 2)


    def test_get_formula(self):

        bc_form_1 = core.BcForm().from_str('abc_a + abc_b')
        self.assertEqual(bc_form_1.get_formula({'abc_a': EmpiricalFormula('C5H10O'), 'abc_b': EmpiricalFormula('C3H5O')}), EmpiricalFormula('C8H15O2'))
        with self.assertRaises(ValueError):
            bc_form_1.get_formula({'abc_a': EmpiricalFormula('C5H10O')})

        bc_form_2 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1+1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        self.assertEqual(bc_form_2.get_formula({'abc_a': EmpiricalFormula('C5H10O'), 'abc_b': EmpiricalFormula('C3H5O')}), EmpiricalFormula('C8H13O'))

        bc_form_3 = core.BcForm().from_str('2 * aa')
        bc_form_3.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(bc_form_3.get_formula(), EmpiricalFormula('C12H26N4O6'))

        bc_form_5 = core.BcForm().from_str('abc_a + abc_b')
        bc_form_5.subunits.append(bc_form_3)
        self.assertEqual(bc_form_5.get_formula({'abc_a': EmpiricalFormula('C5H10O'), 'abc_b': EmpiricalFormula('C3H5O')}), EmpiricalFormula('C20H41N4O8'))

        bc_form_6 = core.BcForm().from_str('2 * aa')
        bc_form_6.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        bc_form_6.subunits.append(bc_form_3)
        self.assertEqual(bc_form_6.get_formula(), EmpiricalFormula('C24H52N8O12'))

        bc_form_7 = core.BcForm().from_str('abc_a + abc_b')
        bc_form_7.set_subunit_attribute('abc_a', 'formula', 'C5H10O')
        with self.assertRaises(ValueError):
            bc_form_7.get_formula()

    def test_get_mol_wt(self):

        bc_form_1 = core.BcForm().from_str('abc_a + abc_b')
        self.assertAlmostEqual(bc_form_1.get_mol_wt({'abc_a': EmpiricalFormula('C5H10O').get_molecular_weight(), 'abc_b': EmpiricalFormula('C3H5O').get_molecular_weight()}), 143, places=0)
        with self.assertRaises(ValueError):
            bc_form_1.get_mol_wt({'abc_a': EmpiricalFormula('C5H10O').get_molecular_weight()})

        bc_form_2 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1+1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        self.assertAlmostEqual(bc_form_2.get_mol_wt({'abc_a': EmpiricalFormula('C5H10O').get_molecular_weight(), 'abc_b': EmpiricalFormula('C3H5O').get_molecular_weight()}), 125, places=0)

        bc_form_3 = core.BcForm().from_str('3 * aa')
        bc_form_3.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        self.assertAlmostEqual(bc_form_3.get_mol_wt(), 483.543, places=3)

        bc_form_5 = core.BcForm().from_str('abc_a + abc_b')
        bc_form_5.subunits.append(bc_form_3)
        self.assertAlmostEqual(bc_form_5.get_mol_wt({'abc_a': EmpiricalFormula('C5H10O').get_molecular_weight(), 'abc_b': EmpiricalFormula('C3H5O').get_molecular_weight()}), EmpiricalFormula('C26H54N6O11').get_molecular_weight())

        bc_form_6 = core.BcForm().from_str('2 * aa')
        bc_form_6.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        bc_form_6.subunits.append(bc_form_3)
        self.assertEqual(bc_form_6.get_mol_wt(), EmpiricalFormula('C30H65N10O15').get_molecular_weight())

        bc_form_7 = core.BcForm().from_str('abc_a + abc_b')
        bc_form_7.set_subunit_attribute('abc_a', 'formula', 'C5H10O')
        with self.assertRaises(ValueError):
            bc_form_7.get_mol_wt()

    def test_get_charge(self):

        bc_form_1 = core.BcForm().from_str('abc_a + abc_b')
        self.assertEqual(bc_form_1.get_charge({'abc_a': 1, 'abc_b': -1}), 0)
        with self.assertRaises(ValueError):
            bc_form_1.get_mol_wt({'abc_a': 1})

        bc_form_2 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1+1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        self.assertEqual(bc_form_2.get_charge({'abc_a': 1, 'abc_b': -1}), -1)

        bc_form_3 = core.BcForm().from_str('2 * aa')
        bc_form_3.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        self.assertEqual(bc_form_3.get_charge(), 2)

        bc_form_5 = core.BcForm().from_str('abc_a + abc_b')
        bc_form_5.subunits.append(bc_form_3)
        self.assertEqual(bc_form_5.get_charge({'abc_a': 1, 'abc_b': -1}), 2)

        bc_form_6 = core.BcForm().from_str('2 * aa')
        bc_form_6.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        bc_form_6.subunits.append(bc_form_3)
        self.assertEqual(bc_form_6.get_charge(), 4)

        bc_form_7 = core.BcForm().from_str('abc_a + abc_b')
        bc_form_7.set_subunit_attribute('abc_a', 'charge', 1)
        with self.assertRaises(ValueError):
            bc_form_7.get_charge()

    def test_validate(self):

        bc_form_1 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        self.assertEqual(len(bc_form_1.validate()), 0)

        bc_form_2 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_c(1)-2O1 | l-displaced-atom: abc_d(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        self.assertEqual(len(bc_form_2.validate()), 2)

        bc_form_3 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a(2)-2O1 | l-displaced-atom: abc_b(2)-2H1 | r-bond-atom: abc_b(3)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        self.assertEqual(len(bc_form_3.validate()), 3)

        bc_form_4 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a-2O1 | l-displaced-atom: abc_a-2H1 | r-bond-atom: abc_b-3C1 | r-displaced-atom: abc_b-3H1 | r-displaced-atom: abc_b-3O1]')
        self.assertEqual(len(bc_form_4.validate()), 0)

        bc_form_5 = core.BcForm().from_str('2 * abc_a + abc_b | x-link: [l-bond-atom: abc_a-2O1 | l-displaced-atom: abc_a-2H1 | r-bond-atom: abc_b-3C1 | r-displaced-atom: abc_b-3H1 | r-displaced-atom: abc_b-3O1]')
        self.assertEqual(len(bc_form_5.validate()), 2)

        bc_form_6 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_c(1)-2O1 | l-displaced-atom: abc_d(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        bc_form_6.subunits.append(bc_form_3)
        self.assertEqual(len(bc_form_6.validate()), 5)

    def test_is_equal(self):

        bc_form_1 = core.BcForm().from_str('abc_a + abc_a + 3 * abc_b')
        bc_form_2 = core.BcForm().from_str('3 * abc_b + 2 * abc_a')
        bc_form_3 = core.BcForm().from_str('abc_a + abc_a + 3 * abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        bc_form_4 = core.BcForm().from_str('abc_a + abc_a + 3 * abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | r-bond-atom: abc_b(1)-3C1 | l-displaced-atom: abc_a(1)-2H1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')
        bc_form_5 = core.BcForm().from_str('abc_a + abc_a + 3 * abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1]')
        bc_form_6 = core.BcForm().from_str('2 * abc_a')
        bc_form_7 = core.BcForm().from_str('4 * abc_a + 3 * abc_b')

        bc_form_8 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a-2O1 | l-displaced-atom: abc_a-2H1 | r-bond-atom: abc_b-3C1 | r-displaced-atom: abc_b-3H1 | r-displaced-atom: abc_b-3O1]')
        bc_form_9 = core.BcForm().from_str('abc_a + abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]')

        bc_form_10 = core.BcForm().from_str('3 * abc_b + 2 * abc_a')
        bc_form_10.subunits.append(bc_form_3)
        bc_form_11 = core.BcForm().from_str('3 * abc_b + 2 * abc_a')
        bc_form_11.subunits.append(bc_form_1)
        bc_form_12 = core.BcForm().from_str('3 * abc_b + 2 * abc_a')
        bc_form_12.subunits.append(core.BcForm().from_str('abc_a + abc_a + 3 * abc_b | x-link: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]'))

        self.assertTrue(bc_form_1.is_equal(bc_form_1))
        self.assertFalse(bc_form_1.is_equal('form'))
        self.assertTrue(bc_form_1.is_equal(bc_form_2))
        self.assertFalse(bc_form_1.is_equal(bc_form_3))
        self.assertTrue(bc_form_3.is_equal(bc_form_4))
        self.assertFalse(bc_form_3.is_equal(bc_form_5))
        self.assertFalse(bc_form_1.is_equal(bc_form_6))
        self.assertFalse(bc_form_1.is_equal(bc_form_7))

        self.assertTrue(bc_form_8.is_equal(bc_form_9))

        self.assertFalse(bc_form_10.is_equal(bc_form_11))
        self.assertTrue(bc_form_10.is_equal(bc_form_12))

    def test_get_subunit_attribute(self):

        bc_form = core.BcForm().from_str('aa')
        bc_form_sub = core.BcForm().from_str('def')
        bc_form.subunits.append(bc_form_sub)
        bc_form.subunits[0].structure = bpforms.ProteinForm().from_str('AA')
        self.assertEqual(bc_form.get_subunit_attribute('aa', 'stoichiometry'), 1)
        self.assertTrue(bc_form.get_subunit_attribute('aa', 'structure').is_equal(bpforms.ProteinForm().from_str('AA')))
        with self.assertRaises(ValueError):
            bc_form.get_subunit_attribute('bb', 'stoichiometry')
        with self.assertRaises(ValueError):
            bc_form.get_subunit_attribute('aa', 'invalidattr')
        self.assertEqual(bc_form.get_subunit_attribute('aa', 'formula'), bpforms.ProteinForm().from_str('AA').get_formula())
        self.assertEqual(bc_form.get_subunit_attribute('aa', 'mol_wt'), bpforms.ProteinForm().from_str('AA').get_mol_wt())
        self.assertEqual(bc_form.get_subunit_attribute('aa', 'charge'), bpforms.ProteinForm().from_str('AA').get_charge())

    def test_set_subunit_attribute(self):

        bc_form = core.BcForm().from_str('aa')
        bc_form_sub = core.BcForm().from_str('def')
        bc_form.subunits.append(bc_form_sub)
        bc_form.set_subunit_attribute('aa', 'structure', bpforms.ProteinForm().from_str('AA'))
        self.assertTrue(bc_form.subunits[0].structure.is_equal(bpforms.ProteinForm().from_str('AA')))
        bc_form.set_subunit_attribute('aa', 'stoichiometry', 2)
        self.assertEqual(bc_form.subunits[0].stoichiometry, 2)
        with self.assertRaises(ValueError):
            bc_form.set_subunit_attribute('bb', 'stoichiometry', 4)
        with self.assertRaises(ValueError):
            bc_form.set_subunit_attribute('aa', 'invalidattr', 'b')

        bc_form.set_subunit_attribute('aa', 'structure', None)
        bc_form.set_subunit_attribute('aa', 'formula', None)
        bc_form.set_subunit_attribute('aa', 'mol_wt', 3.9)
        self.assertTrue(bc_form.subunits[0].mol_wt == 3.9)
        bc_form.set_subunit_attribute('aa', 'charge', 1)
        self.assertTrue(bc_form.subunits[0].charge == 1)
        bc_form.set_subunit_attribute('aa', 'formula', EmpiricalFormula('CH4'))
        self.assertTrue(bc_form.subunits[0].formula == EmpiricalFormula('CH4'))

    def test_get_structure(self):

        # no crosslink
        bc_form_1 = core.BcForm().from_str('2*a')
        self.assertTrue(len(bc_form_1.validate())==0)
        bc_form_1.set_subunit_attribute('a', 'structure', bpforms.ProteinForm().from_str('A'))
        self.assertEqual(OpenBabelUtils.export(bc_form_1.get_structure(), 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)O.C[C@H]([NH3+])C(=O)O')

        # mini "homodimer" AA
        # linking C[C@H]([NH3+])C(=O)[O-] and C[C@H]([NH3+])C(=O)[O-]
        bc_form_2 = core.BcForm().from_str('2*a | x-link: [l-bond-atom: a(1)-1C8 | l-displaced-atom: a(1)-1O10 | l-displaced-atom: a(1)-1H10 | r-bond-atom: a(2)-1N4-1 | r-displaced-atom: a(2)-1H4+1 | r-displaced-atom: a(2)-1H4]')
        self.assertTrue(len(bc_form_2.validate())==0)
        bc_form_2.set_subunit_attribute('a', 'structure', bpforms.ProteinForm().from_str('A'))
        self.assertEqual(OpenBabelUtils.export(bc_form_2.get_structure(), 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)O')

        # mini "heterodimer" AG
        # linking C[C@H]([NH3+])C(=O)[O-] and C([NH3+])C(=O)[O-]
        bc_form_3 = core.BcForm().from_str('a+g | x-link: [l-bond-atom: a-1C8 | l-displaced-atom: a-1O10 | l-displaced-atom: a(1)-1H10 | r-bond-atom: g-1N5-1 | r-displaced-atom: g-1H5+1 | r-displaced-atom: g-1H5]')
        self.assertTrue(len(bc_form_3.validate())==0)
        bc_form_3.set_subunit_attribute('a', 'structure', bpforms.ProteinForm().from_str('A'))
        bc_form_3.set_subunit_attribute('g', 'structure', bpforms.ProteinForm().from_str('G'))
        self.assertEqual(OpenBabelUtils.export(bc_form_3.get_structure(), 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)NCC(=O)O')

        # a more realistic example AGGA, where subunits are composed of multiple monomers
        # linking C[C@H]([NH3+])C(=O)NCC(=O)[O-] and C([NH3+])C(=O)N[C@@H](C)C(=O)[O-]
        bc_form_4 = core.BcForm().from_str('ag+ga | x-link: [l-bond-atom: ag-2C2 | l-displaced-atom: ag-2O1 | l-displaced-atom: ag-2H1 | r-bond-atom: ga-1N5-1 | r-displaced-atom: ga-1H5+1 | r-displaced-atom: ga-1H5]')
        self.assertTrue(len(bc_form_4.validate())==0)
        bc_form_4.set_subunit_attribute('ag', 'structure', bpforms.ProteinForm().from_str('AG'))
        bc_form_4.set_subunit_attribute('ga', 'structure', bpforms.ProteinForm().from_str('GA'))
        self.assertEqual(OpenBabelUtils.export(bc_form_4.get_structure(), 'smiles', options=[]), 'C[C@H]([NH3+])C(=O)NCC(=O)NCC(=O)N[C@@H](C)C(=O)O')

        # a more realistic example ACCMGAGA, where subunits are composed of multiple monomers
        # linking ACCM and 2*GA
        bc_form_5 = core.BcForm().from_str('accm+2*ga | x-link: [l-bond-atom: accm-4C11 | l-displaced-atom: accm-4O13 | l-displaced-atom: accm-4H13 | r-bond-atom: ga(1)-1N5-1 | r-displaced-atom: ga(1)-1H5+1 | r-displaced-atom: ga(1)-1H5] | x-link: [l-bond-atom: ga(1)-2C8 | l-displaced-atom: ga(1)-2O10 | l-displaced-atom: ga(1)-2H10 | r-bond-atom: ga(2)-1N5-1 | r-displaced-atom: ga(2)-1H5+1 | r-displaced-atom: ga(2)-1H5]')
        self.assertTrue(len(bc_form_5.validate())==0)
        bc_form_5.set_subunit_attribute('accm', 'structure', bpforms.ProteinForm().from_str('ACCM'))
        bc_form_5.set_subunit_attribute('ga', 'structure', bpforms.ProteinForm().from_str('GA'))
        self.assertEqual(OpenBabelUtils.export(bc_form_5.get_structure(), 'smiles', options=['canonical']), OpenBabelUtils.export(bpforms.ProteinForm().from_str('ACCMGAGA').get_structure()[0], 'smiles', options=['canonical']))

        # mini "heterodimer" + small molecule AG+CH4
        bc_form_6 = core.BcForm().from_str('a+g+small | x-link: [l-bond-atom: a-1C8 | l-displaced-atom: a-1O10 | l-displaced-atom: a-1H10 | r-bond-atom: g-1N5-1 | r-displaced-atom: g-1H5+1 | r-displaced-atom: g-1H5] | x-link: [l-bond-atom: g-1C4 | l-displaced-atom: g-1H4 | r-bond-atom: small-1C1 | r-displaced-atom: small-1H1 ]')
        self.assertTrue(len(bc_form_6.validate())==0)
        bc_form_6.set_subunit_attribute('a', 'structure', bpforms.ProteinForm().from_str('A'))
        bc_form_6.set_subunit_attribute('g', 'structure', bpforms.ProteinForm().from_str('G'))
        ob_mol = openbabel.OBMol()
        conversion = openbabel.OBConversion()
        conversion.SetInFormat('smi')
        conversion.ReadString(ob_mol, 'C')
        bc_form_6.set_subunit_attribute('small', 'structure', ob_mol)
        self.assertEqual(OpenBabelUtils.export(bc_form_6.get_structure(), 'smiles', options=['canonical']), OpenBabelUtils.export(bpforms.ProteinForm().from_str('AA').get_structure()[0], 'smiles', options=['canonical']))

    def test_export(self):
        bc_form_3 = core.BcForm().from_str('a+g | x-link: [l-bond-atom: a-1C8 | l-displaced-atom: a-1O10 | l-displaced-atom: a(1)-1H10 | r-bond-atom: g-1N5-1 | r-displaced-atom: g-1H5+1 | r-displaced-atom: g-1H5]')
        self.assertTrue(len(bc_form_3.validate())==0)
        bc_form_3.set_subunit_attribute('a', 'structure', bpforms.ProteinForm().from_str('A'))
        bc_form_3.set_subunit_attribute('g', 'structure', bpforms.ProteinForm().from_str('G'))
        self.assertEqual(bc_form_3.export(), 'C[C@H]([NH3+])C(=O)NCC(=O)O')
