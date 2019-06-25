""" bcforms command line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

import cement
import bcforms
import bcforms.core


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "bcforms"
        arguments = [
            (['-v', '--version'], dict(action='version', version=bcforms.__version__)),
        ]

    @cement.ex(hide=True)
    def _default(self):
        self._parser.print_help()


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'bcforms'
        base_controller = 'base'
        handlers = [
            BaseController,
        ]


def main():
    with App() as app:
        app.run()
