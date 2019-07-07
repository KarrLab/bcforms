""" BcForms

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

import itertools
import lark
import pkg_resources
from wc_utils.util.chem import EmpiricalFormula


class Atom(object):
    """ Atom in a crosslink

    Attributes:
        subunit (:obj:`str`): id of subunit
        subunit_idx (:obj:`int`): index of the subunit for homomers
        element (:obj:`str`): code of the element
        position (:obj:`int`): position of the atom within the compound
        monomer (:obj:`int`): index of parent monomer
        charge (:obj:`int`, optional): charge of the atom

    """

    def __init__(self, subunit, element, position, monomer, charge=0, subunit_idx=None):
        """

        Args:
            subunit (:obj:`str`): id of subunit
            element (:obj:`str`): code of the element
            position (:obj:`int`): position of the atom within the compound
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
            subunits (:obj:`list` of :obj:`Subunit`, optional): subunit composition of the complex
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
            :obj:`list` of :obj:`Subunit`: subunits

        """
        return self._subunits

    @subunits.setter
    def subunits(self, value):
        """ Set the subunits

        Args:
            value (:obj:`list` of :obj:`Subunit`): subunits

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
            s += str(subunit['stoichiometry']) + ' * ' + subunit['id'] + ' + '
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
                return [x for x in args if type(x) == dict]

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

            self.subunits.append(new_subunit)

        self.clean()

        return self

    def clean(self):
        """ Clean up the subunits and the crosslinks

        For example, convert `1 * a + 1 * a` to `2 * a`

        """
        subunits_cleaned = []
        subunit_unique_ids = []
        for subunit in self.subunits:
            id = subunit['id']
            if id not in subunit_unique_ids:
                subunit_unique_ids.append(id)
                subunits_cleaned.append(subunit)
            else:
                for subunit_cleaned in subunits_cleaned:
                    if subunit_cleaned['id'] == id:
                        subunit_cleaned['stoichiometry'] += subunit['stoichiometry']
                        break

        self.subunits = subunits_cleaned

    def get_formula(self, subunit_formulas):
        """ Get the empirical formula

        Args:
            subunit_formulas (:obj:`dict`): dictionary of subunit ids and empirical formulas

        Returns:
            :obj:`EmpiricalFormula`: the empirical formula of the BcForm

        Raises:
            :obj:`ValueError`: subunit formulas does not include all subunits

        """

        formula = EmpiricalFormula()

        # subunits
        for subunit in self.subunits:
            if subunit['id'] not in subunit_formulas:
                raise ValueError('subunit_formulas must include all subunits')
            else:
                formula += subunit_formulas[subunit['id']] * subunit['stoichiometry']
        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.left_displaced_atoms, crosslink.right_displaced_atoms):
                formula[atom.element] -= 1
        return formula

    def get_mol_wt(self, subunit_mol_wts):
        """ Get the molecular weight

        Args:
            subunit_formulas (:obj:`dict`): dictionary of subunit ids and molecular weights

        Returns:
            :obj:`float`: the molecular weight of the BcForm

        Raises:
            :obj:`ValueError`: subunit_mol_wts does not include all subunits

        """
        mol_wt = 0.0
        # subunits
        for subunit in self.subunits:
            if subunit['id'] not in subunit_mol_wts:
                raise ValueError('subunit_mol_wts must include all subunits')
            else:
                mol_wt += subunit_mol_wts[subunit['id']] * subunit['stoichiometry']
        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.left_displaced_atoms, crosslink.right_displaced_atoms):
                mol_wt -= EmpiricalFormula(atom.element).get_molecular_weight()

        return mol_wt

    def get_charge(self, subunit_charges):
        """ Get the total charge

        Args:
            subunit_formulas (:obj:`dict`): dictionary of subunit ids and charges

        Returns:
            :obj:`int`: the total charge of the BcForm

        Raises:
            :obj:`ValueError`: subunit_charges does not include all subunits

        """
        charge = 0

        # subunits
        for subunit in self.subunits:
            if subunit['id'] not in subunit_charges:
                raise ValueError('subunit_charges must include all subunits')
            else:
                charge += subunit_charges[subunit['id']] * subunit['stoichiometry']

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
        atom_types = ['left_bond_atoms', 'left_displaced_atoms', 'right_bond_atoms', 'right_displaced_atoms']
        for i_crosslink, crosslink in enumerate(self.crosslinks):
            for atom_type in atom_types:
                for i_atom, atom in enumerate(getattr(crosslink, atom_type)):
                    # check if subunit is present
                    if atom.subunit not in [subunit['id'] for subunit in self.subunits]:
                        errors.append("'{}[{}]' of crosslink {} must belong to a subunit in self.subunits".format(
                            atom_type, i_atom, i_crosslink + 1))
                    # check subunit index
                    elif atom.subunit_idx is None:
                        if next(subunit for subunit in self.subunits if subunit['id'] == atom.subunit)['stoichiometry'] > 1:
                            errors.append("crosslink {} contains multiple subunit '{}', so the subunit_idx of atom '{}[{}]' cannot be None".format(
                            i_crosslink + 1, atom.subunit, atom_type, i_atom))
                    elif atom.subunit_idx > next(subunit for subunit in self.subunits if subunit['id'] == atom.subunit)['stoichiometry']:
                        errors.append("'{}[{}]' of crosslink {} must belong to a subunit whose index is "
                                      "valid in terms of the stoichiometry of the subunit".format(
                                          atom_type, i_atom, i_crosslink + 1))

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
        if sorted(sorted(d.items()) for d in self.subunits) != sorted(sorted(d.items()) for d in other.subunits):
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
