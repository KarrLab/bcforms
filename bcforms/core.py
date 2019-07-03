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

        """


        formula = EmpiricalFormula()

        # subunits
        for subunit in self.subunits:
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
        """
        mol_wt = 0.0
        # subunits
        for subunit in self.subunits:
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
        """
        charge = 0
        # subunits
        for subunit in self.subunits:
            charge += subunit_charges[subunit['id']] * subunit['stoichiometry']
        # crosslinks
        for crosslink in self.crosslinks:
            for atom in itertools.chain(crosslink.left_displaced_atoms, crosslink.right_displaced_atoms):
                charge -= atom.charge

        return charge
