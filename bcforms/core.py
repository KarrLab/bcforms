""" BcForms

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

import itertools
import lark
import openbabel
import pkg_resources
from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
import bpforms
import bpforms.core
import bpforms.alphabet.protein

class Subunit(object):
    """ Subunit in a BcForm macromolecular complex

    Attributes:
        id (:obj:`str`): id of the subunit
        stoichiometry (:obj:`int`): stoichiometry of the subunit
        structure (:obj:`bpforms.BpForm` or :obj:`openbabel.OBMol`, optional): structure of the subunit
        formula (:obj:`EmpiricalFormula`, optional): formula of the subunit
        mol_wt (:obj:`float`, optional): molecular weight of the subunit
        charge (:obj:`int`, optional): charge of the subunit
    """

    def __init__(self, id, stoichiometry, structure=None, formula=None, mol_wt=None, charge=None):
        """

        Args:
            id (:obj:`str`): id of the subunit
            stoichiometry (:obj:`int`): stoichiometry of the subunit
            structure (:obj:`bpforms.BpForm` or :obj:`openbabel.OBMol`, optional): structure of the subunit
            formula (:obj:`EmpiricalFormula`, optional): formula of the subunit
            mol_wt (:obj:`float`, optional): molecular weight of the subunit
            charge (:obj:`int`, optional): charge of the subunit
        """
        self.id = id
        self.stoichiometry = stoichiometry
        self.structure = structure
        if structure is None:
            self.formula = formula
            self.charge = charge
            if formula is None:
                self.mol_wt = mol_wt

    @property
    def id(self):
        """ Get the id of the subunit

        Returns:
            :obj:`str`: id of the subunit

        """
        return self._id

    @id.setter
    def id(self, value):
        """ Set the id of the subunit

        Args:
            value (:obj:`str`): id of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`str`
        """
        if not isinstance(value, str):
            raise ValueError('`value` must be an instance of `str`')
        self._id = value

    @property
    def stoichiometry(self):
        """ Get the stoichiometry of the subunit

        Returns:
            :obj:`int`: stoichiometry of the subunit

        """
        return self._stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, value):
        """ Set the stoichiometry of the subunit

        Args:
            value (:obj:`int`): stoichiometry of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`int`
        """
        if not isinstance(value, int):
            raise ValueError('`value` must be an instance of `int`')
        self._stoichiometry = value

    @property
    def structure(self):
        """ Get the structure of the subunit

        Returns:
            :obj:`bpforms.BpForm` or :obj:`openbabel.OBMol` or None: structure of the subunit

        """
        return self._structure

    @structure.setter
    def structure(self, value):
        """ Set the structure of the subunit

        * setting structure will automaticall set formula, mol_wt, charge

        Args:
            value (:obj:`bpforms.BpForm` or :obj:`openbabel.OBMol` or :obj:`str` (SMILES-encoded string) or None): structure of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`bpforms.BpForm` or :obj:`openbabel.OBMol` or None
            :obj:`ValueError`: unable to convert smiles-encoded string to structure
        """
        if not isinstance(value, bpforms.BpForm) and not isinstance(value, openbabel.OBMol) and not isinstance(value, str) and value is not None:
            raise ValueError('`value` must be an instance of `bpforms.BpForm` or `openbabel.OBMol` or None')

        if isinstance(value, str):
            ob_mol = openbabel.OBMol()
            conversion = openbabel.OBConversion()
            conversion.SetInFormat('smi')
            if not conversion.ReadString(ob_mol, value):
                raise ValueError('unable to convert smiles-encoded string to structure')
            self._structure = ob_mol
        else:
            self._structure = value

        if isinstance(self._structure, openbabel.OBMol):
            self._formula = OpenBabelUtils.get_formula(self._structure)
            self._mol_wt = self.formula.get_molecular_weight()
            self._charge = self._structure.GetTotalCharge()
        elif isinstance(self._structure, bpforms.BpForm):
            self._formula = self._structure.get_formula()
            self._mol_wt = self._structure.get_mol_wt()
            self._charge = self._structure.get_charge()

    @property
    def formula(self):
        """ Get the empirical formula of the subunit

        Returns:
            :obj:`EmpiricalFormula` or None: formula of the subunit

        """
        return self._formula

    @formula.setter
    def formula(self, value):
        """ Set the formula of the subunit

        Args:
            value (:obj:`EmpiricalFormula` or :obj:`str` (string representation of the formula) None): formula of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`EmpiricalFormula` or None
            :obj:`ValueError`: if formula already set by setting structure attribute
        """
        if not isinstance(value, EmpiricalFormula) and not isinstance(value, str) and value is not None:
            raise ValueError(':obj:`value` is not an instance of :obj:`EmpiricalFormula` or :obj:`str` or None')

        if self.structure is not None:
            raise ValueError('formula already set by setting structure attribute')

        if isinstance(value, str):
            self._formula = EmpiricalFormula(value)
        else:
            self._formula = value

        if isinstance(self._formula, EmpiricalFormula):
            self._mol_wt = self._formula.get_molecular_weight()


    @property
    def mol_wt(self):
        """ Get the molecular weight of the subunit

        Returns:
            :obj:`float` or None: molecular weight of the subunit

        """
        return self._mol_wt

    @mol_wt.setter
    def mol_wt(self, value):
        """ Set the molecular weight of the subunit

        Args:
            value (:obj:`float` or :obj:`int` or :obj:`None`): molecular weight of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`float` or :obj:`int` or None
            :obj:`ValueError`: if mol_wt already set by setting structure attribute or formula attribute
            :obj:`ValueError`: if mol_wt is not non-negative
        """
        if not isinstance(value, float) and not isinstance(value, int) and value is not None:
            raise ValueError(':obj:`value` is not an instance of :obj:`float` or :obj:`int` or None')

        if self.formula is not None:
            raise ValueError('mol_wt already set by setting structure attribute or formula attribute')

        if isinstance(value, int):
            value = float(value)

        if isinstance(value, float):
            if value < 0:
                raise ValueError('mol_wt must be non-negative')

        self._mol_wt = value

    @property
    def charge(self):
        """ Get the charge of the subunit

        Returns:
            :obj:`int` or None: charge of the subunit

        """
        return self._charge

    @charge.setter
    def charge(self, value):
        """ Set the charge of the subunit

        Args:
            value (:obj:`int` or :obj:`None`): charge of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`int` or None
            :obj:`ValueError`: if charge already set by setting structure attribute
        """
        if not isinstance(value, int) and value is not None:
            raise ValueError(':obj:`value` is not an instance of :obj:`int` or None')

        if self.structure is not None:
            raise ValueError('charge already set by setting structure attribute')

        self._charge = value

    def __str__(self):
        return str(self.stoichiometry) + ' * ' + self.id

    def is_equal(self, other):
        """ Check if two Subunits are semantically equal

        * Check id and stoichiometry; do not check structure yet

        Args:
            other (:obj:`Subunit`): another Subunit

        Returns:
            :obj:`bool`: :obj:`True`, if the Subunits are semantically equal

        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        attrs = ['id', 'stoichiometry']

        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        return True

    def get_formula(self, formula=None):
        """ Get the empirical formula

        Args:
            formula (:obj:`EmpiricalFormula` or :obj:`None`): Subunit empirical formula per copy

        Returns:
            :obj:`EmpiricalFormula` or None: the empirical formula of the Subunit

        """

        if formula is not None:
            self.formula = formula

        if self.formula is not None:
            return self.formula * self.stoichiometry

        return None

    def get_mol_wt(self, mol_wt=None):
        """ Get the molecular weight

        Args:
            mol_wt (:obj:`float` or :obj:`None`): Subunit molecular weight per copy

        Returns:
            :obj:`float` or None: the molecular weight of the Subunit
        """

        if mol_wt is not None:
            self.mol_wt = mol_wt

        if self.mol_wt is not None:
            return self.mol_wt * self.stoichiometry

        return None

    def get_charge(self, charge=None):
        """ Get the total charge

        Args:
            charge (:obj:`int` or :obj:`None`): Subunit charge per copy

        Returns:
            :obj:`int` or None: the total charge of the Subunit
        """

        if charge is not None:
            self.charge = charge

        if self.charge is not None:
            return self.charge * self.stoichiometry

        return None

    def get_structure(self):
        """ Get an Open Babel molecule of the structure

        Returns:
            :obj:`tuple`:
                * :obj:`openbabel.OBMol`: Open Babel molecule of the structure
                * :obj:`dict` of obj:`dict`: dictionary which maps :obj:`subunit_idx` to
                    atom_maps

        Raises:
            :obj:`ValueError`: Subunit structure is :obj:`None`
        """

        if self.structure is None:
            raise ValueError('Structure is None')

        # join the subunits
        mol = openbabel.OBMol()

        subunit_atom_map = {}
        subunit_idx = 1
        for i in range(self.stoichiometry):

            # get structure
            atom_map = {}
            if isinstance(self.structure, openbabel.OBMol):
                structure = self.structure
                atom_map[1] = {}
                atom_map[1]['monomer'] = {}
                for i_atom in range(structure.NumAtoms()):
                    atom_map[1]['monomer'][i_atom+1] = i_atom+1
            else:
                # structure is a BpForm object
                structure, atom_map = self.structure.get_structure()

            num_atoms = structure.NumAtoms()
            total_atoms = sum(sum(len(y) for y in x.values()) for x in atom_map.values())
            # print(num_atoms, total_atoms)

            mol += structure
            for monomer in atom_map.values():
                for atom_type in monomer.values():
                    for i_atom, atom in atom_type.items():
                        atom_type[i_atom] = atom + num_atoms*(subunit_idx-1)

            subunit_atom_map[subunit_idx] = atom_map
            subunit_idx += 1

        return mol, subunit_atom_map

    def export(self, format='smiles', options=[]):
        """ Export the structure to string

        Args:
            format (:obj:`str`, optional): export format
            options (:obj:`list`, optional): export options

        Returns:
            :obj:`str`: exported string representation of the structure

        """
        if self.structure is None:
            return ''

        return OpenBabelUtils.export(self.get_structure()[0], format=format, options=options)

class Atom(object):
    """ Atom in a crosslink

    Attributes:
        subunit (:obj:`str`): id of subunit
        subunit_idx (:obj:`int`): index of the subunit for homomers
        element (:obj:`str`): code of the element
        position (:obj:`int`): SMILES position of the atom within the compound
        monomer (:obj:`int`): index of parent monomer
        charge (:obj:`int`): charge of the atom
        component_type (:obj:`str`): type of component the atom belongs to:
            either 'monomer' or 'backbone'

    """

    def __init__(self, subunit, element, position, monomer, charge=0, subunit_idx=None, component_type=None):
        """

        Args:
            subunit (:obj:`str`): id of subunit
            element (:obj:`str`): code of the element
            position (:obj:`int`): SMILES position of the atom within the compound
            monomer (:obj:`int`): index of parent monomer
            charge (:obj:`int`, optional): charge of the atom
            subunit_idx (:obj:`int`, optional): index of the subunit for homomers
            component_type (:obj:`str`, optional): type of component the atom belongs to:
                either 'monomer' or 'backbone'
        """

        self.subunit = subunit
        self.subunit_idx = subunit_idx
        self.element = element
        self.position = position
        self.monomer = monomer
        self.charge = charge
        if component_type == 'm':
            self.component_type = 'monomer'
        elif component_type == 'b':
            self.component_type = 'backbone'
        elif component_type is None:
            self.component_type = 'monomer'
        else:
            self.component_type = component_type

    @property
    def subunit(self):
        """ Get the subunit that the atom belongs to

        Returns:
            :obj:`str`: subunit

        """
        return self._subunit

    @subunit.setter
    def subunit(self, value):
        """ Set the subunit that the atom belongs to

        Args:
            value (:obj:`str`): subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`str`
        """
        if not isinstance(value, str):
            raise ValueError('`value` must be an instance of `str`')
        self._subunit = value

    @property
    def subunit_idx(self):
        """ Get the index of the homomer of the subunit that the atom belongs to

        Returns:
            :obj:`int`: subunit_idx or None

        """
        return self._subunit_idx

    @subunit_idx.setter
    def subunit_idx(self, value):
        """ Set the index of the homomer of the subunit that the atom belongs to

        Args:
            value (:obj:`int`): subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not None or a positive integer
        """
        if value is not None and (not isinstance(value, int) or value < 1):
            raise ValueError('`value` must be a None or a positive integer')
        self._subunit_idx = value

    @property
    def element(self):
        """ Get the element of the atom

        Returns:
            :obj:`str`: element

        """
        return self._element

    @element.setter
    def element(self, value):
        """ Set the element of the atom

        Args:
            value (:obj:`str`): element

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`str`
        """
        if not isinstance(value, str):
            raise ValueError('`value` must be an instance of `str`')
        self._element = value

    @property
    def position(self):
        """ Get the position of the atom in the compound

        Returns:
            :obj:`int`: position

        """
        return self._position

    @position.setter
    def position(self, value):
        """ Set the position of the atom in the compound

        Args:
            value (:obj:`int`): position

        Raises:
            :obj:`ValueError`: if :obj:`value` is not a positive :obj:`int`
        """
        if not isinstance(value, int) or value < 1:
            raise ValueError('`value` must be a positive integer')
        self._position = value

    @property
    def monomer(self):
        """ Get the position in the subunit of the monomer that the atom belongs to

        Returns:
            :obj:`int`: monomer position

        """
        return self._monomer

    @monomer.setter
    def monomer(self, value):
        """ Set the position in the subunit of the monomer that the atom belongs to

        Args:
            value (:obj:`int`): monomer position

        Raises:
            :obj:`ValueError`: if `value` is not a positive integer
        """
        if not isinstance(value, int) or value < 1:
            raise ValueError('`value` must be a positive integer')
        self._monomer = value

    @property
    def charge(self):
        """ Get the charge of the atom

        Returns:
            :obj:`int`: charge

        """
        return self._charge

    @charge.setter
    def charge(self, value):
        """ Set the charge of the atom

        Args:
            value (:obj:`int`): charge

        Raises:
            :obj:`ValueError`: if `value` is not an integer
        """
        if not isinstance(value, int):
            raise ValueError('`value` must be an integer')
        self._charge = value

    @property
    def component_type(self):
        """ Get the type of component the atom belongs to

        Returns:
            :obj:`str`: component type

        """
        return self._component_type

    @component_type.setter
    def component_type(self, value):
        """ Set the type of component the atom belongs to

        Raises:
            :obj:`ValueError`: component_type must be either 'monomer' or 'backbone'

        """
        if value not in ['monomer', 'backbone']:
            raise ValueError('`component_type` must be either "monomer" or "backbone"')
        else:
            self._component_type = value

    def __str__(self):
        """ Generate a string representation

        Returns:
            :obj:`str`: string representation
        """

        if self.charge == 0:
            charge = ''
        else:
            charge = '{:+d}'.format(self.charge)

        if self.subunit_idx is None:
            subunit_idx = ''
        else:
            subunit_idx = '(' + str(self.subunit_idx) + ')'
        return '{}{}-{}{}{}{}'.format(self.subunit, subunit_idx, self.monomer, self.element, self.position, charge)

    def is_equal(self, other):
        """ Check if two atoms are semantically equal (belong to the same subunit/monomer and
        have the same element, position, and charge)

        Args:
            other (:obj:`Atom`): another atom

        Returns:
            :obj:`bool`: :obj:`True`, if the atoms are semantically equal

        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        attrs = ['subunit', 'element', 'position', 'monomer', 'charge']

        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        self_subunit_idx = self.subunit_idx if self.subunit_idx is not None else 1
        other_subunit_idx = other.subunit_idx if other.subunit_idx is not None else 1
        if self_subunit_idx != other_subunit_idx:
            return False

        return True


class Crosslink(object):
    """ A crosslink between subunits

    Attributes:
        l_bond_atoms (:obj:`list` of :obj:`Atom`): atoms from the left subunit that bond with the right subunit
        r_bond_atoms (:obj:`list` of :obj:`Atom`): atoms from the right subunit that bond with the left subunit
        l_displaced_atoms (:obj:`list` of :obj:`Atom`): atoms from the left subunit displaced by the crosslink
        r_displaced_atoms (:obj:`list` of :obj:`Atom`): atoms from the right subunit displaced by the crosslink
    """

    def __init__(self, l_bond_atoms=None, r_bond_atoms=None, l_displaced_atoms=None, r_displaced_atoms=None):
        """

        Args:
            l_bond_atoms (:obj:`list`): atoms from the left subunit that bond with the right subunit
            r_bond_atoms (:obj:`list`): atoms from the right subunit that bond with the left subunit
            l_displaced_atoms (:obj:`list`): atoms from the left subunit displaced by the crosslink
            r_displaced_atoms (:obj:`list`): atoms from the right subunit displaced by the crosslink

        """
        if l_bond_atoms is None:
            self.l_bond_atoms = []
        else:
            self.l_bond_atoms = l_bond_atoms

        if r_bond_atoms is None:
            self.r_bond_atoms = []
        else:
            self.r_bond_atoms = r_bond_atoms

        if l_displaced_atoms is None:
            self.l_displaced_atoms = []
        else:
            self.l_displaced_atoms = l_displaced_atoms

        if r_bond_atoms is None:
            self.r_displaced_atoms = []
        else:
            self.r_displaced_atoms = r_bond_atoms

    @property
    def l_bond_atoms(self):
        """ Get the left bond atoms

        Returns:
            :obj:`list` of :obj:`Atom`: left bond atoms

        """
        return self._l_bond_atoms

    @l_bond_atoms.setter
    def l_bond_atoms(self, value):
        """ Set the left bond atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): left bond atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._l_bond_atoms = value

    @property
    def r_bond_atoms(self):
        """ Get the right bond atoms

        Returns:
            :obj:`list` of :obj:`Atom`: right bond atoms

        """
        return self._r_bond_atoms

    @r_bond_atoms.setter
    def r_bond_atoms(self, value):
        """ Set the right bond atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): right bond atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._r_bond_atoms = value

    @property
    def l_displaced_atoms(self):
        """ Get the left displaced atoms

        Returns:
            :obj:`list` of :obj:`Atom`: left displaced atoms

        """
        return self._l_displaced_atoms

    @l_displaced_atoms.setter
    def l_displaced_atoms(self, value):
        """ Set the left displaced atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): left displaced atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._l_displaced_atoms = value

    @property
    def r_displaced_atoms(self):
        """ Get the right displaced atoms

        Returns:
            :obj:`list` of :obj:`Atom`: right displaced atoms

        """
        return self._r_displaced_atoms

    @r_displaced_atoms.setter
    def r_displaced_atoms(self, value):
        """ Set the right displaced atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): right displaced atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._r_displaced_atoms = value

    def __str__(self):
        """Generate a string representation

        Returns:
            :obj:`str`: string representation
        """
        s = 'x-link: ['
        atom_types = ['l_bond_atoms', 'l_displaced_atoms', 'r_bond_atoms', 'r_displaced_atoms']
        for atom_type in atom_types:
            for atom in getattr(self, atom_type):
                s += ' {}: {} |'.format(atom_type[:-1].replace('_', '-'), str(atom))

        s = s[:-1]+']'
        return s

    def is_equal(self, other):
        """ Check if two crosslinks are semantically equal (have the same bond atoms)

        Args:
            other (:obj:`Crosslink`): another crosslink

        Returns:
            :obj:`bool`: :obj:`True`, if the crosslinks are semantically equal

        """

        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        attrs = ['l_bond_atoms', 'l_displaced_atoms', 'r_bond_atoms', 'r_displaced_atoms']

        for attr in attrs:
            self_atoms = getattr(self, attr)
            other_atoms = getattr(other, attr)
            if len(self_atoms) != len(other_atoms):
                return False
            for self_atom, other_atom in zip(self_atoms, other_atoms):
                if not self_atom.is_equal(other_atom):
                    return False

        return True


class BcForm(object):
    """ A form of a macromolecular complex

    Attributes:
        subunits (:obj:`list` of :obj:`Subunit`): subunit composition of the complex
        crosslinks (:obj:`list` :obj:`Crosslink`): crosslinks in the complex

    """

    def __init__(self, subunits=None, crosslinks=None):
        """

        Args:
            subunits (:obj:`list` of :obj:`Subunit` or :obj:`BcForm`, optional): subunit composition of the complex
            crosslinks (:obj:`list` of :obj:`Crosslink`, optional): crosslinks in the complex

            _parser (:obj:`lark.Lark`): lark grammar parser used in `from_str`
        """
        if subunits is None:
            self.subunits = []
        else:
            self.subunits = subunits

        if crosslinks is None:
            self.crosslinks = []
        else:
            self.crosslinks = crosslinks

    @property
    def subunits(self):
        """ Get the subunits

        Returns:
            :obj:`list` of :obj:`Subunit` or :obj:`BcForm`: subunits

        """
        return self._subunits

    @subunits.setter
    def subunits(self, value):
        """ Set the subunits

        Args:
            value (:obj:`list` of :obj:`Subunit` or :obj`BcForm`): subunits

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._subunits = value

    @property
    def crosslinks(self):
        """ Get the crosslinks

        Returns:
            :obj:`list` of :obj:`Crosslink`: crosslinks

        """
        return self._crosslinks

    @crosslinks.setter
    def crosslinks(self, value):
        """ Set the crosslinks

        Args:
            value (:obj:`list` of :obj:`Crosslink`): crosslinks

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._crosslinks = value

    def __str__(self):
        """ Generate a string representation

        Returns:
            :obj:`str`: string representation of complex
        """
        s = ''

        # subunits
        for subunit in self.subunits:
            s += str(subunit) + ' + '
        s = s[:-3]

        # crosslinks
        for crosslink in self.crosslinks:
            s += ' | ' + str(crosslink)

        # return string representation
        return s

    # read the grammar file
    # _grammar_filename = 'grammar.lark'
    _grammar_filename = pkg_resources.resource_filename('bcforms', 'grammar.lark')

    with open(_grammar_filename, 'r') as file:
        _parser = lark.Lark(file.read())

    def from_str(self, string):
        """ Set a complex from a string representation

        Args:
            string (:obj:`str`): string representation of a complex

        Returns:
            :obj:`BcForm`: structured BcForm representation of the string
        """

        class ParseTreeTransformer(lark.Transformer):
            # Class that processes the parsetree

            def __init__(self, bc_form):
                super(ParseTreeTransformer, self).__init__()
                self.bc_form = bc_form

            @lark.v_args(inline=True)
            def start(self, *args):
                self.bc_form.subunits = args[0]
                self.bc_form.crosslinks = []
                if len(args) > 2:
                    # exists global attr (crosslink)
                    self.bc_form.crosslinks = list(args[2::2])
                return self.bc_form

            # complex
            @lark.v_args(inline=True)
            def complex(self, *args):
                return [Subunit(id=x['id'], stoichiometry=x['stoichiometry']) for x in args if type(x) == dict]

            @lark.v_args(inline=True)
            def component(self, *args):
                component_dict = {}
                if len(args) < 2:
                    # handle the case where no explicit coefficient
                    component_dict['stoichiometry'] = 1
                    component_dict[args[0][0]] = args[0][1]
                else:
                    # handle the case where optional coefficient is explicitly put
                    component_dict[args[0][0]] = args[0][1]
                    component_dict[args[1][0]] = args[1][1]

                return component_dict

            @lark.v_args(inline=True)
            def coefficient(self, *args):
                return ('stoichiometry', int(args[0].value))

            @lark.v_args(inline=True)
            def subunit(self, *args):
                return ('id', args[0].value)

            # crosslinks
            @lark.v_args(inline=True)
            def global_attr(self, *args):
                return args[0]

            @lark.v_args(inline=True)
            def crosslink(self, *args):
                bond = Crosslink()
                for arg in args:
                    if isinstance(arg, tuple):
                        atom_type, atom = arg
                        atom_type_list = getattr(bond, atom_type+"s")
                        atom_type_list.append(atom)
                return bond

            @lark.v_args(inline=True)
            def crosslink_atom(self, *args):
                num_optional_args = 0
                atom_type = args[0][1]
                subunit = args[2][1]
                if args[3][0] == 'subunit_idx':
                    subunit_idx = int(args[3][1])
                else:
                    subunit_idx = None
                    num_optional_args += 1
                monomer = int(args[4-num_optional_args][1])
                element = args[5-num_optional_args][1]
                position = int(args[6-num_optional_args][1])
                if len(args) > 7-num_optional_args:
                    if args[7-num_optional_args][0] == 'atom_component_type':
                        atom_component_type = args[7-num_optional_args][1]
                    else:
                        atom_component_type = None
                        num_optional_args += 1
                else:
                    atom_component_type = None
                if len(args) > 8-num_optional_args:
                    if args[8-num_optional_args][0] == 'atom_charge':
                        charge = int(args[8-num_optional_args][1])
                    else:
                        charge = 0
                else:
                    charge = 0

                return (atom_type, Atom(subunit=subunit, subunit_idx=subunit_idx, element=element,
                                        position=position, monomer=monomer, charge=charge, component_type=atom_component_type))

            @lark.v_args(inline=True)
            def crosslink_atom_type(self, *args):
                return ('crosslink_atom_type', args[0].value + '_' + args[1].value + '_atom')

            @lark.v_args(inline=True)
            def monomer_position(self, *args):
                return ('monomer_position', int(args[0].value))

            @lark.v_args(inline=True)
            def subunit_idx(self, *args):
                return ('subunit_idx', int(args[0].value[1:-1]))

            @lark.v_args(inline=True)
            def atom_element(self, *args):
                return ('atom_element', args[0].value)

            @lark.v_args(inline=True)
            def atom_position(self, *args):
                return ('atom_position', int(args[0].value))

            @lark.v_args(inline=True)
            def atom_charge(self, *args):
                return ('atom_charge', args[0].value)

            @lark.v_args(inline=True)
            def atom_component_type(self, *args):
                return ('atom_component_type', args[0].value)

        tree = self._parser.parse(string)
        # print(tree.pretty())
        parse_tree_transformer = ParseTreeTransformer(self)
        bc_form = parse_tree_transformer.transform(tree)
        bc_form.clean()
        return bc_form

    def from_set(self, subunits):
        """ Set the subunits from a list of subunits

        Note: this method does not support crosslinks

        Args:
            subunits: (:obj:`list`): list representation of a complex. For example::

                [
                    {'id': 'ABC_A', 'stoichiometry': 2},
                    {'id': 'ABC_B', 'stoichiometry': 3},
                ]

        Returns:
            :obj:`BcForm`: this complex

        Raises:
            :obj:`ValueError`: subunit has no 'id' key
            :obj:`ValueError`: subunit has no 'stoichiometry' key
        """
        self.subunits = []
        self.crosslinks = []

        for subunit in subunits:
            new_subunit = {}

            # process id of subunit
            if 'id' in subunit:
                new_subunit['id'] = subunit['id']
            else:
                raise ValueError('`subunit` has no `id`')

            # process stoichiometry of subunit
            if 'stoichiometry' in subunit:
                new_subunit['stoichiometry'] = subunit['stoichiometry']
            else:
                raise ValueError('`subunit` has no `stoichiometry`')

            self.subunits.append(Subunit(id=new_subunit['id'], stoichiometry=new_subunit['stoichiometry']))

        self.clean()

        return self

    def clean(self):
        """ Clean up the subunits and the crosslinks

        For example, convert `1 * a + 1 * a` to `2 * a`

        """
        subunits_cleaned = []
        subunit_unique_ids = []
        for subunit in self.subunits:
            if isinstance(subunit, Subunit):
                id = subunit.id
                if id not in subunit_unique_ids:
                    subunit_unique_ids.append(id)
                    subunits_cleaned.append(subunit)
                else:
                    next(subunit_cleaned for subunit_cleaned in subunits_cleaned if subunit_cleaned.id == id).stoichiometry += subunit.stoichiometry
            elif isinstance(subunit, BcForm):
                subunit.clean()
                subunits_cleaned.append(subunit)

        self.subunits = subunits_cleaned

    def get_formula(self, subunit_formulas=None):
        """ Get the empirical formula

        * If user wants to calculate formula of nested BcForm, where some subunits
          are BcForm objects, then the subunit BcForms must be able to calculate
          its own formula through structure

        Args:
            subunit_formulas (:obj:`dict` or :obj:`None`): dictionary of subunit ids and empirical formulas

        Returns:
            :obj:`EmpiricalFormula`: the empirical formula of the BcForm

        Raises:
            :obj:`ValueError`: subunit formulas does not include all subunits
            :obj:`ValueError`: Not all subunits have defined formula
        """

        formula = EmpiricalFormula()

        # subunits
        if subunit_formulas is None:
            for subunit in self.subunits:
                if subunit.get_formula() is None:
                    raise ValueError('Not all subunits have defined formula')
                formula += subunit.get_formula()
        else:
            for subunit in self.subunits:
                if isinstance(subunit, BcForm):
                    formula += subunit.get_formula()
                else:
                    if subunit.id not in subunit_formulas:
                        raise ValueError('subunit_formulas must include all subunits')
                    else:
                        formula += subunit.get_formula(formula=subunit_formulas[subunit.id])

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.l_displaced_atoms, crosslink.r_displaced_atoms):
                formula[atom.element] -= 1
        return formula

    def get_mol_wt(self, subunit_mol_wts=None):
        """ Get the molecular weight

        * If user wants to calculate molecular weight of nested BcForm, where
          some subunits are BcForm objects, then the subunit BcForms must be able
          to calculate its own molecular weight through structure

        Args:
            subunit_formulas (:obj:`dict` or :obj:`None`): dictionary of subunit ids and molecular weights

        Returns:
            :obj:`float`: the molecular weight of the BcForm

        Raises:
            :obj:`ValueError`: subunit_mol_wts does not include all subunits
            :obj:`ValueError`: Not all subunits have defined molecular weight
        """
        mol_wt = 0.0

        # subunits
        if subunit_mol_wts is None:
            for subunit in self.subunits:
                if subunit.get_mol_wt() is None:
                    raise ValueError('Not all subunits have defined molecular weight')
                mol_wt += subunit.get_mol_wt()
        else:
            for subunit in self.subunits:
                if isinstance(subunit, BcForm):
                    mol_wt += subunit.get_mol_wt()
                else:
                    if subunit.id not in subunit_mol_wts:
                        raise ValueError('subunit_mol_wts must include all subunits')
                    else:
                        mol_wt += subunit.get_mol_wt(mol_wt=subunit_mol_wts[subunit.id])

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.l_displaced_atoms, crosslink.r_displaced_atoms):
                mol_wt -= EmpiricalFormula(atom.element).get_molecular_weight()

        return mol_wt

    def get_charge(self, subunit_charges=None):
        """ Get the total charge

        * If user wants to calculate charge of nested BcForm, where
          some subunits are BcForm objects, then the subunit BcForms must be able
          to calculate its own charge through structure

        Args:
            subunit_formulas (:obj:`dict` or :obj:`None`): dictionary of subunit ids and charges

        Returns:
            :obj:`int`: the total charge of the BcForm

        Raises:
            :obj:`ValueError`: subunit_charges does not include all subunits
            :obj:`ValueError`: Not all subunits have defined charge
        """
        charge = 0

        # subunits
        if subunit_charges is None:
            for subunit in self.subunits:
                if subunit.get_charge() is None:
                    raise ValueError('Not all subunits have defined charge')
                charge += subunit.get_charge()
        else:
            for subunit in self.subunits:
                if isinstance(subunit, BcForm):
                    charge += subunit.get_charge()
                else:
                    if subunit.id not in subunit_charges:
                        raise ValueError('subunit_charges must include all subunits')
                    else:
                        charge += subunit.get_charge(charge=subunit_charges[subunit.id])

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.l_displaced_atoms, crosslink.r_displaced_atoms):
                charge -= atom.charge

        # return the total charge
        return charge

    def validate(self):
        """ Check if the BcForm is valid

        * Check if the crosslinking subunit is in the subunit list and if the `subunit_idx` is valid

        Returns:
            :obj:`list` of :obj:`str`: list of errors, if any

        """
        errors = []

        # crosslinks
        self_subunits_subunits = [subunit for subunit in self.subunits if isinstance(subunit,Subunit)]
        self_subunits_bcforms = [subunit for subunit in self.subunits if isinstance(subunit,BcForm)]

        atom_types = ['l_bond_atoms', 'l_displaced_atoms', 'r_bond_atoms', 'r_displaced_atoms']
        for i_crosslink, crosslink in enumerate(self.crosslinks):
            for atom_type in atom_types:
                for i_atom, atom in enumerate(getattr(crosslink, atom_type)):
                    # check if subunit is present
                    if atom.subunit not in [subunit.id for subunit in self_subunits_subunits]:
                        errors.append("'{}[{}]' of crosslink {} must belong to a subunit in self.subunits".format(
                            atom_type, i_atom, i_crosslink + 1))
                    # check subunit index
                    elif atom.subunit_idx is None:
                        if next(subunit for subunit in self_subunits_subunits if subunit.id == atom.subunit).stoichiometry > 1:
                            errors.append("crosslink {} contains multiple subunit '{}', so the subunit_idx of atom '{}[{}]' cannot be None".format(
                            i_crosslink + 1, atom.subunit, atom_type, i_atom))
                    elif atom.subunit_idx > next(subunit for subunit in self_subunits_subunits if subunit.id == atom.subunit).stoichiometry:
                        errors.append("'{}[{}]' of crosslink {} must belong to a subunit whose index is "
                                      "valid in terms of the stoichiometry of the subunit".format(
                                          atom_type, i_atom, i_crosslink + 1))

        for self_subunits_bcform in self_subunits_bcforms:
            errors.extend(self_subunits_bcform.validate())

        return errors

    def is_equal(self, other):
        """ Check if two complexes are semantically equal (same subunits and crosslinks)

        Args:
            other (:obj:`BcForm`): another complex

        Returns:
            :obj:`bool`: :obj:`True`, if the complexes are semantically equal

        """

        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        # test subunits
        if len(self.subunits) != len(other.subunits):
            return False
        for self_subunit in self.subunits:
            found = False
            for other_subunit in other.subunits:
                if self_subunit.is_equal(other_subunit):
                    found = True
                    break
            if not found:
                return False

        # test crosslinks
        if len(self.crosslinks) != len(other.crosslinks):
            return False
        for self_crosslink in self.crosslinks:
            found = False
            for other_crosslink in other.crosslinks:
                if self_crosslink.is_equal(other_crosslink):
                    found = True
                    break
            if not found:
                return False

        return True

    def get_subunit_attribute(self, subunit_id, attribute):
        """ Set attribute (stoichiometry, structure) of subunit by id

        Args:
            subunit_id (:obj:`str`): id of subunit
            attribute (:obj:`str`): attribute to set

        Returns:
            :obj:`int` for stoichiometry, :obj:`bpforms.BpForm`, :obj:`openbabel.OBMol`, or None for structure

        Raises:
            :obj:`ValueError`: No Subunit with subunit_id
            :obj:`ValueError`: Invalid attribute
        """

        subunit = next((subunit for subunit in self.subunits if isinstance(subunit, Subunit) and subunit.id == subunit_id), None)
        if subunit is None:
            raise ValueError('No Subunit with subunit_id')

        if attribute not in ['stoichiometry', 'structure', 'formula', 'mol_wt', 'charge']:
            raise ValueError('Invalid attribute')

        return getattr(subunit, attribute)

    def set_subunit_attribute(self, subunit_id, attribute, value):
        """ Set attribute (stoichiometry, structure) of subunit by id

        Args:
            subunit_id (:obj:`str`): id of subunit
            attribute (:obj:`str`): attribute to set
            value (:obj:`int` for stoichiometry, :obj:`bpforms.BpForm`, :obj:`openbabel.OBMol`, or None for structure): value

        Raises:
            :obj:`ValueError`: No Subunit with subunit_id
            :obj:`ValueError`: Invalid attribute
        """

        subunit = next((subunit for subunit in self.subunits if isinstance(subunit, Subunit) and subunit.id == subunit_id), None)
        if subunit is None:
            raise ValueError('No Subunit with subunit_id')

        if attribute not in ['stoichiometry', 'structure', 'formula', 'mol_wt', 'charge']:
            raise ValueError('Invalid attribute')

        setattr(subunit, attribute, value)

    def get_structure(self):
        """ Get an Open Babel molecule of the structure

        Returns:
            :obj:`openbabel.OBMol`: Open Babel molecule of the structure
        """
        mol = openbabel.OBMol()

        atom_maps = []
        n_atoms = [0]

        # subunits
        for i_subunit, subunit in enumerate(self.subunits):
            structure, atom_map = subunit.get_structure()
            mol += structure

            n_atoms.append(n_atoms[-1]+structure.NumAtoms())

            for subunit_map in atom_map.values():
                for monomer in subunit_map.values():
                    for atom_type in monomer.values():
                        for i_atom, atom in atom_type.items():
                            atom_type[i_atom] = atom+n_atoms[i_subunit]
            atom_maps.append(atom_map)

        # mol.AddHydrogens()

        # print(atom_maps)
        # for i in range(mol.NumAtoms()):
        #     print(mol.GetAtom(i+1), mol.GetAtom(i+1).GetAtomicNum())

        bonding_hydrogens = []
        # crosslinks
        # get the atoms
        crosslinks_atoms = []
        for crosslink in self.crosslinks:
            crosslink_atoms = {}
            crosslinks_atoms.append(crosslink_atoms)
            for atom_type in ['l_bond_atoms', 'r_bond_atoms', 'l_displaced_atoms', 'r_displaced_atoms']:
                crosslink_atoms[atom_type] = []
                for atom_md in getattr(crosslink, atom_type):
                    i_subunit = [i for i in range(len(self.subunits)) if self.subunits[i].id == atom_md.subunit][0]
                    subunit_idx = 1 if atom_md.subunit_idx is None else atom_md.subunit_idx
                    atom = mol.GetAtom(atom_maps[i_subunit][subunit_idx][atom_md.monomer][atom_md.component_type][atom_md.position])
                    if atom_md.element == 'H' and atom.GetAtomicNum() != 1:
                        atom = get_hydrogen_atom(atom, bonding_hydrogens, (i_subunit, subunit_idx-1, atom_md.monomer-1, atom_md.component_type))
                    crosslink_atoms[atom_type].append((atom, atom_md.charge))

        # print(OpenBabelUtils.export(mol, format='smiles', options=[]))

        # make the crosslink bonds
        for atoms in crosslinks_atoms:

            for atom, atom_charge in itertools.chain(atoms['l_displaced_atoms'], atoms['r_displaced_atoms']):
                if atom:
                    assert mol.DeleteAtom(atom, True)

            for (l_atom, l_atom_charge), (r_atom, r_atom_charge) in zip(atoms['l_bond_atoms'], atoms['r_bond_atoms']):
                bond = openbabel.OBBond()
                bond.SetBegin(l_atom)
                bond.SetEnd(r_atom)
                bond.SetBondOrder(1)
                assert mol.AddBond(bond)

                if l_atom_charge:
                    l_atom.SetFormalCharge(l_atom.GetFormalCharge() + l_atom_charge)

                if r_atom_charge:
                    r_atom.SetFormalCharge(r_atom.GetFormalCharge() + r_atom_charge)

        return mol

    def export(self, format='smiles', options=[]):
        """ Export the structure to string

        Args:
            format (:obj:`str`, optional): export format
            options (:obj:`list`, optional): export options

        Returns:
            :obj:`str`: exported string representation of the structure

        """
        return OpenBabelUtils.export(self.get_structure(), format=format, options=options)


def get_hydrogen_atom(parent_atom, bonding_hydrogens, i_monomer):
    """ Get a hydrogen atom attached to a parent atom
    Args:
        parent_atom (:obj:`openbabel.OBAtom`): parent atom
        bonding_hydrogens (:obj:`list`): hydrogens that have already been gotten
        i_monomer (:obj:`int`): index of parent monomer in sequence
    Returns:
        :obj:`openbabel.OBAtom`: hydrogen atom
    """
    for other_atom in openbabel.OBAtomAtomIter(parent_atom):
        if other_atom.GetAtomicNum() == 1:
            tmp = (i_monomer, other_atom.GetIdx())
            if tmp not in bonding_hydrogens:  # hydrogen
                bonding_hydrogens.append(tmp)
                return other_atom
    return None
