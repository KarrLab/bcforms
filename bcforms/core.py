""" bcforms

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

import pkg_resources
import lark
from wc_utils.util.chem import EmpiricalFormula
import itertools

class Atom(object):
    """ atom in crosslink

    Attributes:
        subunit (:obj:`str`): id of subunit
        subunit_idx (:obj:`int`): index of the subunit for homo..mers
        element (:obj:`str`): code of the element
        position (:obj:`int`): position of the atom within the compound
        monomer (:obj:`int`): index of parent monomer
        charge (:obj:`int`, optional): charge of the atom

    """

    def __init__(self, subunit, subunit_idx, element, position, monomer, charge=0):
        """

        Args:
            subunit (:obj:`str`): id of subunit
            subunit_idx (:obj:`int`, optional): index of the subunit for homo..mers
            element (:obj:`str`): code of the element
            position (:obj:`int`): position of the atom within the compound
            monomer (:obj:`int`): index of parent monomer
            charge (:obj:`int`, optional): charge of the atom

        """

        self.subunit = subunit
        self.subunit_idx = subunit_idx
        self.element = element
        self.position = position
        self.monomer = monomer
        self.charge = charge

    @property
    def subunit(self):
        """ get the subunit the atom belongs to

        Returns:
            :obj:`str`: subunit

        """
        return self._subunit

    @subunit.setter
    def subunit(self, value):
        """ set the subunit the atom belongs to

        Args:
            :obj:`str`: subunit

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `str`
        """
        if not isinstance(value, str):
            raise ValueError('`value` must be an instance of `str`')
        self._subunit = value

    @property
    def subunit_idx(self):
        """ get the subunit_idx of the subunit the atom belongs to

        Returns:
            :obj:`int`: subunit_idx

        """
        return self._subunit_idx

    @subunit_idx.setter
    def subunit_idx(self, value):
        """ set the subunit the atom belongs to

        Args:
            :obj:`int`: subunit

        Raises:
            :obj:`ValueError`: if `value` is not a positive int'
        """
        if (not isinstance(value, int) or value<1):
            raise ValueError('`value` must be a positive int')
        self._subunit_idx = value

    @property
    def element(self):
        """ get the element of the atom

        Returns:
            :obj:`str`: element

        """
        return self._element

    @element.setter
    def element(self, value):
        """ set the element of the atom

        Args:
            :obj:`str`: element

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `str`
        """
        if not isinstance(value, str):
            raise ValueError('`value` must be an instance of `str`')
        self._element = value

    @property
    def position(self):
        """ get the position of the atom in the compound

        Returns:
            :obj:`int`: position

        """
        return self._position

    @position.setter
    def position(self, value):
        """ set the position of the atom in the compound

        Args:
            :obj:`int`: position

        Raises:
            :obj:`ValueError`: if `value` is not a positive int'
        """
        if (not isinstance(value, int) or value<1):
            raise ValueError('`value` must be a positive int')
        self._position = value

    @property
    def monomer(self):
        """ get the position in the subunit of the monomer the atom belongs to

        Returns:
            :obj:`int`: monomer position

        """
        return self._monomer

    @monomer.setter
    def monomer(self, value):
        """ set the position in the subunit of the monomer the atom belongs to

        Args:
            :obj:`int`: monomer position

        Raises:
            :obj:`ValueError`: if `value` is not a positive int'
        """
        if (not isinstance(value, int) or value<1):
            raise ValueError('`value` must be a positive int')
        self._monomer = value

    @property
    def charge(self):
        """ get the charge of the atom

        Returns:
            :obj:`int`: charge

        """
        return self._charge

    @charge.setter
    def charge(self, value):
        """ set the charge of the atom

        Args:
            :obj:`int`: charge

        Raises:
            :obj:`ValueError`: if `value` is not an int'
        """
        if not isinstance(value, int):
            raise ValueError('`value` must be an int')
        self._charge = value

    def __str__(self):

        if self.charge == 0:
            charge = ''
        else:
            charge = '%+d' % self.charge

        return '{}({})-{}{}{}{}'.format(self.subunit, self.subunit_idx, self.monomer, self.element, self.position, charge)

    def is_equal(self, other):
        """ check if two Atoms are semantically equal

        Args:
            other (:obj:`Atom`): another Atom

        Returns:
            :obj:`bool`: :obj:`True`, if the Atoms have the same structure

        """
        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        attrs = ['subunit', 'subunit_idx', 'element', 'position', 'monomer', 'charge']

        for attr in attrs:
            if getattr(self, attr) != getattr(other, attr):
                return False

        return True


class Crosslink(object):
    """ crosslink between subunits

    Attributes:
        left_bond_atoms (:obj:`list`): atoms from the left subunit that bond with the right subunit
        right_bond_atoms (:obj:`list`): atoms from the right subunit that bond with the left subunit
        left_displaced_atoms (:obj:`list`): atoms from the left subunit displaced by the crosslink
        right_displaced_atoms (:obj:`list`): atoms from the right subunit displaced by the crosslink
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
        """ get the left_bond_atoms

        Returns:
            :obj:`list`: left_bond_atoms

        """
        return self._left_bond_atoms

    @left_bond_atoms.setter
    def left_bond_atoms(self, value):
        """ set the left_bond_atoms

        Args:
            :obj:`list`: left_bond_atoms

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._left_bond_atoms = value

    @property
    def right_bond_atoms(self):
        """ get the right_bond_atoms

        Returns:
            :obj:`list`: right_bond_atoms

        """
        return self._right_bond_atoms

    @right_bond_atoms.setter
    def right_bond_atoms(self, value):
        """ set the right_bond_atoms

        Args:
            :obj:`list`: right_bond_atoms

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._right_bond_atoms = value

    @property
    def left_displaced_atoms(self):
        """ get the left_displaced_atoms

        Returns:
            :obj:`list`: left_displaced_atoms

        """
        return self._left_displaced_atoms

    @left_displaced_atoms.setter
    def left_displaced_atoms(self, value):
        """ set the left_displaced_atoms

        Args:
            :obj:`list`: left_displaced_atoms

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._left_displaced_atoms = value

    @property
    def right_displaced_atoms(self):
        """ get the right_displaced_atoms

        Returns:
            :obj:`list`: right_displaced_atoms

        """
        return self._right_displaced_atoms

    @right_displaced_atoms.setter
    def right_displaced_atoms(self, value):
        """ set the right_displaced_atoms

        Args:
            :obj:`list`: right_displaced_atoms

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._right_displaced_atoms = value

    def __str__(self):
        s = 'crosslink: ['
        atom_types = ['left_bond_atoms', 'left_displaced_atoms', 'right_bond_atoms', 'right_displaced_atoms']
        for atom_type in atom_types:
            for atom in getattr(self, atom_type):
                s += ' {}: {} |'.format('-'.join(atom_type.split('_'))[:-1], str(atom))

        s = s[:-1]+']'
        return s

    def is_equal(self, other):
        """ check if two Crosslinks are semantically equal

        Args:
            other (:obj:`Crosslink`): another Crosslink

        Returns:
            :obj:`bool`: :obj:`True`, if the Crosslinks have the same structure

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
            for i in range(len(self_atoms)):
                if not self_atoms[i].is_equal(other_atoms[i]):
                    return False

        return True

class BcForm(object):
    """ Biocomplex form

    Attributes:
        subunits (:obj:`list`): subunits composition of the Biocomplex
        crosslinks (:obj:`list`): crosslinks in the Biocomplex

    """

    def __init__(self, subunits=None, crosslinks=None):
        """

        Args:
            subunits (:obj:`list`): subunits composition of the Biocomplex
            crosslinks (:obj:`list`): crosslinks in the Biocomplex

            _parser (:obj:`lark.Lark`): lark grammar parser used in from_str
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
        """ get the subunits of BcForm

        Returns:
            :obj:`list`: subunits

        """
        return self._subunits

    @subunits.setter
    def subunits(self, value):
        """ set the subunits of BcForm

        Args:
            :obj:`list`: subunits

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._subunits = value

    @property
    def crosslinks(self):
        """ get the crosslinks of BcForm

        Returns:
            :obj:`list`: crosslinks

        """
        return self._crosslinks

    @crosslinks.setter
    def crosslinks(self, value):
        """ set the crosslinks of BcForm

        Args:
            :obj:`list`: crosslinks

        Raises:
            :obj:`ValueError`: if `value` is not an instance of `list`

        """
        if not isinstance(value, list):
            raise ValueError('`value` must be an instance of `list`')
        self._crosslinks = value

    def __str__(self):
        s = ''

        # subunits
        for i in range(len(self.subunits)-1):
            s += str(self.subunits[i]['stoichiometry']) + ' * '+ self.subunits[i]['id'] + ' + '
        s += str(self.subunits[-1]['stoichiometry']) + ' * '+ self.subunits[-1]['id']

        # crosslinks
        for crosslink in self.crosslinks:
            s += ' | ' + str(crosslink)

        return s

    # read the grammar file
    # _grammar_filename = 'grammar.lark'
    _grammar_filename = pkg_resources.resource_filename('bcforms', 'grammar.lark')

    with open(_grammar_filename, 'r') as file:
        _parser = lark.Lark(file.read())

    def from_str(self, string):
        """ create a Biocomplex from a string representation

        Args:
            string: (:obj:`str`): string representation of a Biocomplex

        Returns:
            :obj:`BcForm`: the BcForm object of the string
        """

        # class that process the parsetree
        class ParseTreeTransformer(lark.Transformer):

            def __init__(self, bc_form):
                super(ParseTreeTransformer, self).__init__()
                self.bc_form = bc_form

            @lark.v_args(inline=True)
            def start(self, *args):
                self.bc_form.subunits = args[0]
                self.bc_form.crosslinks = []
                if len(args)>2:
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
                atom_type = args[0][1]
                subunit = args[2][1]
                subunit_idx = int(args[3][1])
                monomer = int(args[4][1])
                element = args[5][1]
                position = int(args[6][1])
                if len(args)>7:
                    charge = int(args[7][1])
                else:
                    charge = 0

                return (atom_type, Atom(subunit=subunit, subunit_idx=subunit_idx, element=element, position=position, monomer=monomer, charge=charge))


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

    def from_set(self, set):
        """ create a Biocomplex from a set representation. Cannot represent crosslinks

        Args:
            set: (:obj:`list`): set representation of a Biocomplex (ex. [{'id': 'ABC_A', 'stoichiometry': 2}, {'id': 'ABC_B', 'stoichiometry': 3}])

        Returns:
            :obj:`BcForm`: the BcForm object of the string

        Raises:
            :obj:`ValueError`: subunit has no 'id'
            :obj:`ValueError`: subunit has no 'stoichiometry'
        """
        bc_form = BcForm()
        for subunit in set:
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

            bc_form.subunits.append(new_subunit)

        bc_form.clean()

        return bc_form

    def clean(self):
        """ clean up the subunits and the crosslinks

        """
        # convert 1*a + 1*a to 2*a
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

        self.subunits = subunits_cleaned

    def get_formula(self, subunit_formulas):
        """ get the Empirical Formula of the BcForm complex

        Args:
            subunit_formulas (:obj:`dict`): dictionary of subunit ids and empirical formulas

        Returns:
            :obj:`EmpiricalFormula`: the empirical formula of the BcForm

        Raises:
            :obj:`ValueError`: subunit_formulas must include all subunits

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
        """ get the molecular weight of the BcForm complex

        Args:
            subunit_formulas (:obj:`dict`): dictionary of subunit ids and molecular weights

        Returns:
            :obj:`float`: the molecular weight of the BcForm

        Raises:
            :obj:`ValueError`: subunit_mol_wts must include all subunits

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
        """ get the total charges of the BcForm complex

        Args:
            subunit_formulas (:obj:`dict`): dictionary of subunit ids and charges

        Returns:
            :obj:`int`: the total charge of the BcForm

        Raises:
            :obj:`ValueError`: subunit_charges must include all subunits

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

        return charge

    def validate(self):
        """ check if the BcForm is valid
            for now, check if the crosslinking subunit is in the subunit list
            and if the subunit_idx is valid

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
                        errors.append("'{}[{}]' of crosslink {} must belong to a subunit in self.subunits".format(atom_type, i_atom, i_crosslink+1))
                    # check subunit index
                    elif atom.subunit_idx > next(subunit for subunit in self.subunits if subunit['id']==atom.subunit)['stoichiometry']:
                        errors.append("'{}[{}]' of crosslink {} must belong to a subunit whose index is valid in terms of the stoichiometry of the subunit".format(atom_type, i_atom, i_crosslink+1))

        return errors

    def is_equal(self, other):
        """ check if two BcForms are semantically equal

        Args:
            other (:obj:`BcForm`): another BcForm

        Returns:
            :obj:`bool`: :obj:`True`, if the BcForms have the same structure

        """

        if self is other:
            return True
        if self.__class__ != other.__class__:
            return False

        # test subunit
        if len(self.subunits) != len(other.subunits):
            return False
        if sorted(sorted(d.items()) for d in self.subunits) != sorted(sorted(d.items()) for d in other.subunits):
            return False

        # test crosslink
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


# if __name__ == '__main__':
    # atom_1 = Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
    # atom_2 = Atom(subunit='abc', subunit_idx=1, element='H', position=1, monomer=10, charge=0)
    # atom_3 = Atom(subunit='def', subunit_idx=1, element='H', position=1, monomer=10, charge=0)

    # print(atom_1.is_equal(atom_1))
    # print(atom_1.is_equal('atom'))
    # print(atom_1.is_equal(atom_2))
    # print(atom_1.is_equal(atom_3))

    # crosslink_1 = Crosslink()
    # crosslink_1.left_bond_atoms.append(atom_1)
    # crosslink_1.right_bond_atoms.append(atom_2)
    #
    # crosslink_2 = Crosslink()
    # crosslink_2.left_bond_atoms.append(atom_1)
    # crosslink_2.right_bond_atoms.append(atom_2)
    #
    # crosslink_3 = Crosslink()
    # crosslink_3.left_bond_atoms.append(atom_1)
    # crosslink_3.right_bond_atoms.append(atom_2)
    # crosslink_3.right_bond_atoms.append(atom_3)
    #
    # crosslink_4 = Crosslink()
    # crosslink_4.left_bond_atoms.append(atom_1)
    # crosslink_4.right_bond_atoms.append(atom_3)
    #
    #
    # print(crosslink_1.is_equal(crosslink_1))
    # print(crosslink_1.is_equal(atom_1))
    # print(crosslink_1.is_equal(crosslink_2))
    # print(crosslink_1.is_equal(crosslink_3))
    # print(crosslink_1.is_equal(crosslink_4))

    # bc_form_1 = BcForm().from_str('abc_a + abc_a + 3 * abc_b')
    # bc_form_2 = BcForm().from_str('3 * abc_b + 2 * abc_a')
    # bc_form_3 = BcForm().from_str('abc_a + abc_a + 3 * abc_b | crosslink: [left-bond-atom: abc_a(1)-2O1 | left-displaced-atom: abc_a(1)-2H1 | right-bond-atom: abc_b(1)-3C1 | right-displaced-atom: abc_b(1)-3H1 | right-displaced-atom: abc_b(1)-3O1]')
    # bc_form_4 = BcForm().from_str('abc_a + abc_a + 3 * abc_b | crosslink: [left-bond-atom: abc_a(1)-2O1 | right-bond-atom: abc_b(1)-3C1 | left-displaced-atom: abc_a(1)-2H1 | right-displaced-atom: abc_b(1)-3H1 | right-displaced-atom: abc_b(1)-3O1]')
    # bc_form_5 = BcForm().from_str('abc_a + abc_a + 3 * abc_b | crosslink: [left-bond-atom: abc_a(1)-2O1 | left-displaced-atom: abc_a(1)-2H1 | right-bond-atom: abc_b(1)-3C1 | right-displaced-atom: abc_b(1)-3H1]')
    #
    #
    # print(bc_form_1.is_equal(bc_form_1))
    # print(bc_form_1.is_equal('form'))
    # print(bc_form_1.is_equal(bc_form_2))
    # print(bc_form_1.is_equal(bc_form_3))
    # print(bc_form_3.is_equal(bc_form_4))
    # print(bc_form_3.is_equal(bc_form_5))
