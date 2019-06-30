""" bcforms

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

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


    def from_str(self, string):
        """ create a Biocomplex from a string representation

        Args:
            string: (:obj:`str`): string representation of a Biocomplex

        Returns:
            :obj:`BcForm`: the BcForm object of the string
        """
        pass
