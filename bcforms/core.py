""" bcforms

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

import pkg_resources
import lark

class BcForm(object):
    """ Biocomplex form

    Attributes:
        subunits (:obj:`list`): subunits composition of the Biocomplex
        crosslinks (:obj:`list`): crosslinks in the Biocomplex

    """

    def __init__(self, subunits=[], crosslinks=[]):
        """

        Args:
            subunits (:obj:`list`): subunits composition of the Biocomplex
            crosslinks (:obj:`list`): crosslinks in the Biocomplex

            _parser (:obj:`lark.Lark`): lark grammar parser used in from_str
        """
        self.subunits = subunits
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
                    component_dict["stoichiometry"] = 1
                    component_dict[args[0][0]] = args[0][1]
                else:
                    # handle the case where optional coefficient is explicitly put
                    component_dict[args[0][0]] = args[0][1]
                    component_dict[args[1][0]] = args[1][1]

                return component_dict

            @lark.v_args(inline=True)
            def coefficient(self, *args):
                return ("stoichiometry", int(args[0].value))

            @lark.v_args(inline=True)
            def subunit(self, *args):
                return ("id", args[0].value)

            # crosslinks
            @lark.v_args(inline=True)
            def global_attr(self, *args):
                return args[0]

            @lark.v_args(inline=True)
            def crosslink(self, *args):
                return [x for x in args if type(x) == dict]

            @lark.v_args(inline=True)
            def crosslink_atom(self, *args):
                crosslink_dict = {}
                for arg in args:
                    if isinstance(arg, tuple):
                        arg_name, arg_val = arg
                        crosslink_dict[arg_name] = arg_val
                return crosslink_dict


            @lark.v_args(inline=True)
            def crosslink_atom_type(self, *args):
                return ('crosslink_atom_type', args[0].value+"_"+args[1].value)

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
        parse_tree_transformer = ParseTreeTransformer(self)
        return parse_tree_transformer.transform(tree)
