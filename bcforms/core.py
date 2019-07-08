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
    """

    def __init__(self, id, stoichiometry, structure=None):
        """

        Args:
            id (:obj:`str`): id of the subunit
            stoichiometry (:obj:`int`): stoichiometry of the subunit
            structure (:obj:`bpforms.BpForm` or :obj:`openbabel.OBMol`, optional): structure of the subunit

        """
        self.id = id
        self.stoichiometry = stoichiometry
        self.structure = structure

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

        Args:
            value (:obj:`bpforms.BpForm` or :obj:`openbabel.OBMol` or None): structure of the subunit

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`bpforms.BpForm` or :obj:`openbabel.OBMol` or None
        """
        if not isinstance(value, bpforms.BpForm) and not isinstance(value, openbabel.OBMol) and value is not None:
            raise ValueError('`value` must be an instance of `bpforms.BpForm` or `openbabel.OBMol` or None')
        self._structure = value

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

    def get_formula(self, formula=None, from_structure=False):
        """ Get the empirical formula

        Args:
            formula (:obj:`EmpiricalFormula` or None): Subunit empirical formula per copy
            from_structure (:obj:`bool`): True if calculate formula from subunit structure

        Returns:
            :obj:`EmpiricalFormula`: the empirical formula of the Subunit

        Raises:
            :obj:`ValueError`: Attempting to get formula from structure but structure is None

        """

        if formula is None or from_structure is True:
            # get formula from structure
            if self.structure is None:
                raise ValueError('Attempting to get formula from structure but structure is None')
            elif isinstance(self.structure, openbabel.OBMol):
                return OpenBabelUtils.get_formula(self.structure) * self.stoichiometry
            else:
                # structure is a BpForm object
                return self.structure.get_formula() * self.stoichiometry
        else:
            # get formula from formula
            return formula * self.stoichiometry

    def get_mol_wt(self, mol_wt=None, from_structure=False):
        """ Get the molecular weight

        Args:
            mol_wt (:obj:`float` or None): Subunit molecular weight per copy
            from_structure (:obj:`bool`): True if calculate formula from subunit structure

        Returns:
            :obj:`float`: the molecular weight of the Subunit

        Raises:
            :obj:`ValueError`: Attempting to get molecular weight from structure but structure is None

        """

        if mol_wt is None or from_structure is True:
            # get mol_wt from structure
            if self.structure is None:
                raise ValueError('Attempting to get molecular weight from structure but structure is None')
            elif isinstance(self.structure, openbabel.OBMol):
                return OpenBabelUtils.get_formula(self.structure).get_molecular_weight() * self.stoichiometry
            else:
                # structure is a BpForm object
                return self.structure.get_mol_wt() * self.stoichiometry
        else:
            # get mol_wt from formula
            return mol_wt * self.stoichiometry

    def get_charge(self, charge=None, from_structure=False):
        """ Get the total charge

        Args:
            charge (:obj:`int` or None): Subunit charge per copy
            from_structure (:obj:`bool`): True if calculate formula from subunit structure

        Returns:
            :obj:`int`: the total charge of the Subunit

        Raises:
            :obj:`ValueError`: Attempting to get total charge from structure but structure is None

        """

        if charge is None or from_structure is True:
            # get charge from structure
            if self.structure is None:
                raise ValueError('Attempting to get molecular weight from structure but structure is None')
            elif isinstance(self.structure, openbabel.OBMol):
                return self.structure.GetTotalCharge() * self.stoichiometry
            else:
                # structure is a BpForm object
                return self.structure.get_charge() * self.stoichiometry
        else:
            # get charge from formula
            return charge * self.stoichiometry

    def get_structure(self):
        """ Get an OpenBabel molecule of the structure

        Returns:
            :obj:`openbabel.OBMol`: OpenBabel molecule of the structure

        Raises:
            :obj:`ValueError`: Subunit structure is None

        """
        if self.structure is None:
            raise ValueError('Structure is None')
        elif isinstance(self.structure, openbabel.OBMol):
            structure = self.structure
        else:
            # structure is a BpForm object
            structure = self.structure.get_structure()

        mol = openbabel.OBMol()
        for i in range(self.stoichiometry):
            mol += structure
        return mol


class Atom(object):
    """ Atom in a crosslink

    Attributes:
        subunit (:obj:`str`): id of subunit
        subunit_idx (:obj:`int`): index of the subunit for homomers
        element (:obj:`str`): code of the element
        position (:obj:`int`): SMILES position of the atom within the compound
        monomer (:obj:`int`): index of parent monomer
        charge (:obj:`int`, optional): charge of the atom

    """

    def __init__(self, subunit, element, position, monomer, charge=0, subunit_idx=None):
        """

        Args:
            subunit (:obj:`str`): id of subunit
            element (:obj:`str`): code of the element
            position (:obj:`int`): SMILES position of the atom within the compound
            monomer (:obj:`int`): index of parent monomer
            charge (:obj:`int`, optional): charge of the atom
            subunit_idx (:obj:`int`, optional): index of the subunit for homomers
        """

        self.subunit = subunit
        self.subunit_idx = subunit_idx
        self.element = element
        self.position = position
        self.monomer = monomer
        self.charge = charge

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
        left_bond_atoms (:obj:`list` of :obj:`Atom`): atoms from the left subunit that bond with the right subunit
        right_bond_atoms (:obj:`list` of :obj:`Atom`): atoms from the right subunit that bond with the left subunit
        left_displaced_atoms (:obj:`list` of :obj:`Atom`): atoms from the left subunit displaced by the crosslink
        right_displaced_atoms (:obj:`list` of :obj:`Atom`): atoms from the right subunit displaced by the crosslink
    """

    def __init__(self, left_bond_atoms=None, right_bond_atoms=None, left_displaced_atoms=None, right_displaced_atoms=None):
        """

        Args:
            left_bond_atoms (:obj:`list`): atoms from the left subunit that bond with the right subunit
            right_bond_atoms (:obj:`list`): atoms from the right subunit that bond with the left subunit
            left_displaced_atoms (:obj:`list`): atoms from the left subunit displaced by the crosslink
            right_displaced_atoms (:obj:`list`): atoms from the right subunit displaced by the crosslink

        """
        if left_bond_atoms is None:
            self.left_bond_atoms = []
        else:
            self.left_bond_atoms = left_bond_atoms

        if right_bond_atoms is None:
            self.right_bond_atoms = []
        else:
            self.right_bond_atoms = right_bond_atoms

        if left_displaced_atoms is None:
            self.left_displaced_atoms = []
        else:
            self.left_displaced_atoms = left_displaced_atoms

        if right_bond_atoms is None:
            self.right_displaced_atoms = []
        else:
            self.right_displaced_atoms = right_bond_atoms

    @property
    def left_bond_atoms(self):
        """ Get the left bond atoms

        Returns:
            :obj:`list` of :obj:`Atom`: left bond atoms

        """
        return self._left_bond_atoms

    @left_bond_atoms.setter
    def left_bond_atoms(self, value):
        """ Set the left bond atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): left bond atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._left_bond_atoms = value

    @property
    def right_bond_atoms(self):
        """ Get the right bond atoms

        Returns:
            :obj:`list` of :obj:`Atom`: right bond atoms

        """
        return self._right_bond_atoms

    @right_bond_atoms.setter
    def right_bond_atoms(self, value):
        """ Set the right bond atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): right bond atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._right_bond_atoms = value

    @property
    def left_displaced_atoms(self):
        """ Get the left displaced atoms

        Returns:
            :obj:`list` of :obj:`Atom`: left displaced atoms

        """
        return self._left_displaced_atoms

    @left_displaced_atoms.setter
    def left_displaced_atoms(self, value):
        """ Set the left displaced atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): left displaced atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._left_displaced_atoms = value

    @property
    def right_displaced_atoms(self):
        """ Get the right displaced atoms

        Returns:
            :obj:`list` of :obj:`Atom`: right displaced atoms

        """
        return self._right_displaced_atoms

    @right_displaced_atoms.setter
    def right_displaced_atoms(self, value):
        """ Set the right displaced atoms

        Args:
            value (:obj:`list` of :obj:`Atom`): right displaced atoms

        Raises:
            :obj:`ValueError`: if :obj:`value` is not an instance of :obj:`list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._right_displaced_atoms = value

    def __str__(self):
        """Generate a string representation

        Returns:
            :obj:`str`: string representation
        """
        s = 'crosslink: ['
        atom_types = ['left_bond_atoms', 'left_displaced_atoms', 'right_bond_atoms', 'right_displaced_atoms']
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

        attrs = ['left_bond_atoms', 'left_displaced_atoms', 'right_bond_atoms', 'right_displaced_atoms']

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
                    charge = int(args[7-num_optional_args][1])
                else:
                    charge = 0

                return (atom_type, Atom(subunit=subunit, subunit_idx=subunit_idx, element=element,
                                        position=position, monomer=monomer, charge=charge))

            @lark.v_args(inline=True)
            def crosslink_atom_type(self, *args):
                return ('crosslink_atom_type', args[0].value+'_'+args[1].value+'_atom')

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

    def get_formula(self, subunit_formulas=None, from_structure=False):
        """ Get the empirical formula

        * if user wants to calculate formula of nested BcForm, where some subunits
        are BcForm objects, then the subunit BcForms must be able to calculate
        its own formula through structure

        Args:
            subunit_formulas (:obj:`dict` or None): dictionary of subunit ids and empirical formulas
            from_structure (:obj:`bool`): True if calculate formula from subunit structure

        Returns:
            :obj:`EmpiricalFormula`: the empirical formula of the BcForm

        Raises:
            :obj:`ValueError`: subunit formulas does not include all subunits

        """

        formula = EmpiricalFormula()

        # subunits
        # if no subunit_formulas provided, then try to get formula from structure
        # (even if from_structure is False by default)
        # or if subunit_formulas provided and from_structure is True, then
        # try to get formula from structure
        if subunit_formulas is None or from_structure is True:
            for subunit in self.subunits:
                formula += subunit.get_formula(from_structure=True)
        # if subunit_formulas provided and from_structure is False, then
        # get formula from subunit_formulas
        else:
            for subunit in self.subunits:
                if isinstance(subunit, BcForm):
                    formula += subunit.get_formula(from_structure=True)
                else:
                    if subunit.id not in subunit_formulas:
                        raise ValueError('subunit_formulas must include all subunits')
                    else:
                        formula += subunit.get_formula(formula=subunit_formulas[subunit.id], from_structure=False)

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.left_displaced_atoms, crosslink.right_displaced_atoms):
                formula[atom.element] -= 1
        return formula

    def get_mol_wt(self, subunit_mol_wts=None, from_structure=False):
        """ Get the molecular weight

        * if user wants to calculate molecular weight of nested BcForm, where
        some subunits are BcForm objects, then the subunit BcForms must be able
        to calculate its own molecular weight through structure

        Args:
            subunit_formulas (:obj:`dict` or None): dictionary of subunit ids and molecular weights
            from_structure (:obj:`bool`): True if calculate molecular weight from subunit structure

        Returns:
            :obj:`float`: the molecular weight of the BcForm

        Raises:
            :obj:`ValueError`: subunit_mol_wts does not include all subunits

        """
        mol_wt = 0.0

        # subunits
        # if no subunit_mol_wts provided, then try to get mol_wt from structure
        # (even if from_structure is False by default)
        # or if subunit_mol_wts provided and from_structure is True, then
        # try to get mol_wt from structure
        if subunit_mol_wts is None or from_structure is True:
            for subunit in self.subunits:
                mol_wt += subunit.get_mol_wt(from_structure=True)
        # if subunit_mol_wts provided and from_structure is False, then
        # get mol_wt from subunit_formulas
        else:
            for subunit in self.subunits:
                if isinstance(subunit, BcForm):
                    mol_wt += subunit.get_mol_wt(from_structure=True)
                else:
                    if subunit.id not in subunit_mol_wts:
                        raise ValueError('subunit_mol_wts must include all subunits')
                    else:
                        mol_wt += subunit.get_mol_wt(mol_wt=subunit_mol_wts[subunit.id], from_structure=False)

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.left_displaced_atoms, crosslink.right_displaced_atoms):
                mol_wt -= EmpiricalFormula(atom.element).get_molecular_weight()

        return mol_wt

    def get_charge(self, subunit_charges=None, from_structure=False):
        """ Get the total charge

        * if user wants to calculate charge of nested BcForm, where
        some subunits are BcForm objects, then the subunit BcForms must be able
        to calculate its own charge through structure

        Args:
            subunit_formulas (:obj:`dict` or None): dictionary of subunit ids and charges
            from_structure (:obj:`bool`): True if calculate molecular weight from subunit structure

        Returns:
            :obj:`int`: the total charge of the BcForm

        Raises:
            :obj:`ValueError`: subunit_charges does not include all subunits

        """
        charge = 0

        # subunits
        # if no subunit_charges provided, then try to get charge from structure
        # (even if from_structure is False by default)
        # or if subunit_charges provided and from_structure is True, then
        # try to get charge from structure
        if subunit_charges is None or from_structure is True:
            for subunit in self.subunits:
                charge += subunit.get_charge(from_structure=True)
        # if subunit_charges provided and from_structure is False, then
        # get charge from subunit_formulas
        else:
            for subunit in self.subunits:
                if isinstance(subunit, BcForm):
                    charge += subunit.get_charge(from_structure=True)
                else:
                    if subunit.id not in subunit_charges:
                        raise ValueError('subunit_charges must include all subunits')
                    else:
                        charge += subunit.get_charge(charge=subunit_charges[subunit.id], from_structure=False)

        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.left_displaced_atoms, crosslink.right_displaced_atoms):
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

        atom_types = ['left_bond_atoms', 'left_displaced_atoms', 'right_bond_atoms', 'right_displaced_atoms']
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

        if attribute not in ['stoichiometry', 'structure']:
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

        if attribute not in ['stoichiometry', 'structure']:
            raise ValueError('Invalid attribute')

        setattr(subunit, attribute, value)

    def get_structure(self):
        """ Get an OpenBabel molecule of the structure

        Returns:
            :obj:`openbabel.OBMol`: OpenBabel molecule of the structure

        """
        mol = openbabel.OBMol()

        i_atoms = [0]
        n_atoms = []
        n_monomers_atoms = []

        # subunits
        for subunit in self.subunits:
            structure = subunit.get_structure()
            mol += structure
            i_atoms.append(i_atoms[-1]+structure.NumAtoms())
            n_atoms.append(int(structure.NumAtoms()/subunit.stoichiometry))

            n_monomer_atom = [0]
            if isinstance(subunit.structure, bpforms.BpForm):
                for i in range(len(subunit.structure)):
                    # print(OpenBabelUtils.export(subunit.structure[i].structure, 'smiles', options=[]))
                    n_monomer_atom.append(n_monomer_atom[-1]+subunit.structure[i].structure.NumAtoms())

            n_monomers_atoms.append(n_monomer_atom)

        # crosslinks
        # get the atoms
        crosslinks_atoms = []
        for crosslink in self.crosslinks:
            crosslink_atoms = {}
            crosslinks_atoms.append(crosslink_atoms)
            for atom_type in ['left_bond_atoms', 'right_bond_atoms', 'left_displaced_atoms', 'right_displaced_atoms']:
                crosslink_atoms[atom_type] = []
                for atom_md in getattr(crosslink, atom_type):
                    # calculate the index of the crosslinking atom in the molecule
                    i_subunit = [i for i in range(len(self.subunits)) if self.subunits[i].id == atom_md.subunit][0]
                    subunit_idx = 1 if atom_md.subunit_idx is None else atom_md.subunit_idx
                    # print(atom_md.monomer-1)
                    # print(n_monomers_atoms[i_subunit][atom_md.monomer-1])
                    i_atom = i_atoms[i_subunit] + n_atoms[i_subunit]*(subunit_idx-1) + n_monomers_atoms[i_subunit][atom_md.monomer-1] + atom_md.position
                    # get the atom
                    atom = mol.GetAtom(i_atom)
                    # print(atom_type, atom_md.element, atom.GetAtomicNum())
                    crosslink_atoms[atom_type].append(atom)

        # make the crosslink bonds
        for atoms in crosslinks_atoms:
            for l_atom, r_atom in zip(atoms['left_bond_atoms'], atoms['right_bond_atoms']):
                bond = openbabel.OBBond()
                bond.SetBegin(l_atom)
                bond.SetEnd(r_atom)
                bond.SetBondOrder(1)
                assert mol.AddBond(bond)
            for atom in itertools.chain(atoms['left_displaced_atoms'], atoms['right_displaced_atoms']):
                if atom:
                    assert mol.DeleteAtom(atom, True)

        return mol


if __name__ == '__main__':
    print('Example: no crosslink')
    print('should be: C[C@H]([NH3+])C(=O)[O-].C[C@H]([NH3+])C(=O)[O-]')
    bc_form_1 = BcForm().from_str('2*a')
    assert len(bc_form_1.validate())==0
    bc_form_1.set_subunit_attribute('a', 'structure', bpforms.alphabet.protein.ProteinForm().from_str('A'))
    print('   is    :',OpenBabelUtils.export(bc_form_1.get_structure(), 'smiles', options=[]))
    print()

    print('Example: "mini homodimer AA"')
    print('Linking C[C@H]([NH3+])C(=O)[O-] and C[C@H]([NH3+])C(=O)[O-]')
    print('should be: C[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)[O-]')
    bc_form_2 = BcForm().from_str('2*a | crosslink: [left-bond-atom: a(1)-1C8 | left-displaced-atom: a(1)-1O10 | right-bond-atom: a(2)-1N4 | right-displaced-atom: a(2)-1H5 | right-displaced-atom: a(2)-1H6]')
    assert len(bc_form_2.validate())==0
    bc_form_2.set_subunit_attribute('a', 'structure', bpforms.alphabet.protein.ProteinForm().from_str('A'))
    print('   is    :',OpenBabelUtils.export(bc_form_2.get_structure(), 'smiles', options=[]))
    print()

    print('Example: "mini heterodimer AG"')
    print('Linking C[C@H]([NH3+])C(=O)[O-] and C([NH3+])C(=O)[O-]')
    print('should be: C[C@H]([NH3+])C(=O)NCC(=O)[O-]')
    bc_form_3 = BcForm().from_str('a+g | crosslink: [left-bond-atom: a-1C8 | left-displaced-atom: a-1O10 | right-bond-atom: g-1N2 | right-displaced-atom: g-1H3 | right-displaced-atom: g-1H4]')
    assert len(bc_form_3.validate())==0
    bc_form_3.set_subunit_attribute('a', 'structure', bpforms.alphabet.protein.ProteinForm().from_str('A'))
    bc_form_3.set_subunit_attribute('g', 'structure', bpforms.alphabet.protein.ProteinForm().from_str('G'))
    print('   is    :',OpenBabelUtils.export(bc_form_3.get_structure(), 'smiles', options=[]))
    print()

    print('Example: "a more realistic example AGGA, where subunits are composed of multiple monomers"')
    print('Linking C[C@H]([NH3+])C(=O)NCC(=O)[O-] and C([NH3+])C(=O)N[C@@H](C)C(=O)[O-]')
    print('should be: C[C@H]([NH3+])C(=O)NCC(=O)NCC(=O)N[C@@H](C)C(=O)[O-]')
    bc_form_4 = BcForm().from_str('ag+ga | crosslink: [left-bond-atom: ag-2C6 | left-displaced-atom: ag-2O8 | right-bond-atom: ga-1N2 | right-displaced-atom: ga-1H3 | right-displaced-atom: ga-1H4]')
    assert len(bc_form_4.validate())==0
    bc_form_4.set_subunit_attribute('ag', 'structure', bpforms.alphabet.protein.ProteinForm().from_str('AG'))
    bc_form_4.set_subunit_attribute('ga', 'structure', bpforms.alphabet.protein.ProteinForm().from_str('GA'))
    print('   is    :',OpenBabelUtils.export(bc_form_4.get_structure(), 'smiles', options=[]))
