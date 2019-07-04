""" REST JSON API

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-07-03
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bcforms
import bcforms.core
from wc_utils.util.chem import EmpiricalFormula
import flask
import flask_restplus
import flask_restplus.errors
import flask_restplus.fields

# setup app
app = flask.Flask(__name__)


class PrefixMiddleware(object):
    def __init__(self, app, prefix=''):
        self.app = app
        self.prefix = prefix

    def __call__(self, environ, start_response):
        if environ['PATH_INFO'].startswith(self.prefix):
            environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
            environ['SCRIPT_NAME'] = self.prefix
            return self.app(environ, start_response)
        else:
            start_response('404', [('Content-Type', 'text/plain')])
            return ["This url does not belong to the app.".encode()]


app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix='/api')

api = flask_restplus.Api(app,
                         title='bcforms JSON REST API',
                         description='JSON REST API for calculating properties of biocomplex forms',
                         contact='karr@mssm.edu',
                         version=bcforms.__version__,
                         license='MIT',
                         license_url='https://github.com/KarrLab/bcforms/blob/master/LICENSE',
                         doc='/')

bcform_ns = flask_restplus.Namespace('bcform', description='Calculate properties of biocomplex forms')
api.add_namespace(bcform_ns)

bcforms_model = bcform_ns.model('BcForm', {
    'form': flask_restplus.fields.String(required=True, title='BcForm', description='input biocomplex form', example='2 * abc_a + 3 * abc_b'),
    'subunit_formulas': flask_restplus.fields.String(required=False, title='subunit_formulas', description='formulas of the subunits', example='{abc_a:C5H10O, abc_b:C3H5O}'),
    'subunit_mol_wts': flask_restplus.fields.String(required=False, title='subunit_mol_wts', description='molecular weights of the subunits', example='{abc_a:86, abc_b:57}'),
    'subunit_charges': flask_restplus.fields.String(required=False, title='subunit_charges', description='charges of the subunits', example='{abc_a:+1, abc_b:-1}'),
})


@bcform_ns.route("/")
class Bcform(flask_restplus.Resource):

    @bcform_ns.expect(bcforms_model, validate=True)
    def post(self):
        ret = {}

        args = bcform_ns.payload

        # read arguments
        form = args['form']
        arg_subunit_formulas = args.get('subunit_formulas', None)
        arg_subunit_mol_wts = args.get('subunit_mol_wts', None)
        arg_subunit_charges = args.get('subunit_charges', None)

        # validate form
        try:
            bc_form = bcforms.core.BcForm().from_str(form)
        except Exception as error:
            flask_restplus.abort(400, 'Form is invalid', errors={'form': str(error)})

        errors = bc_form.validate()
        if errors:
            flask_restplus.abort(400, 'Form is invalid', errors={'form': '. '.join(errors)})

        ret['form'] = str(bc_form)

        # if subunit_formulas exists, then calculate biocomplex formula
        if arg_subunit_formulas is not None:

            # parse formula
            subunit_formulas = {}
            try:
                for subunit in arg_subunit_formulas[1:-1].split(','):
                    id, formula = subunit.strip().split(':')
                    subunit_formulas[id.strip()] = EmpiricalFormula(formula.strip())
            except Exception as error:
                flask_restplus.abort(400, 'Cannot parse subunit_formulas', errors={'subunit_formulas': str(error)})

            # calculate BcForm formula
            try:
                formula = bc_form.get_formula(subunit_formulas)
            except Exception as error:
                flask_restplus.abort(400, 'Unable to calculate BcForm formula', errors={'subunit_formulas': str(error)})

            ret['formula'] = str(formula)

        # if subunit_mol_wts exists, then calculate biocomplex molecular weight
        if arg_subunit_mol_wts is not None:

            # parse subunit_mol_wts
            subunit_mol_wts = {}
            try:
                for subunit in arg_subunit_mol_wts[1:-1].split(','):
                    id, mol_wt = subunit.strip().split(':')
                    subunit_mol_wts[id.strip()] = float(mol_wt.strip())
            except Exception as error:
                flask_restplus.abort(400, 'Cannot parse subunit_mol_wts', errors={'subunit_mol_wts': str(error)})

            # calculate BcForm molecular weights
            try:
                mol_wt = bc_form.get_mol_wt(subunit_mol_wts)
            except Exception as error:
                flask_restplus.abort(400, 'Unable to calculate BcForm molecular weight', errors={'subunit_mol_wts': str(error)})

            ret['mol_wt'] = str(mol_wt)

        # if subunit_charges exists, then calculate biocomplex total charges
        if arg_subunit_charges is not None:

            # parse subunit_charges
            subunit_charges = {}
            try:
                for subunit in arg_subunit_charges[1:-1].split(','):
                    id, charge = subunit.strip().split(':')
                    subunit_charges[id.strip()] = int(charge.strip())
            except Exception as error:
                flask_restplus.abort(400, 'Cannot parse subunit_charges', errors={'subunit_charges': str(error)})

            # calculate BcForm charge
            try:
                charge = bc_form.get_charge(subunit_charges)
            except Exception as error:
                flask_restplus.abort(400, 'Unable to calculate BcForm charges', errors={'subunit_charges': str(error)})

            ret['charge'] = str(charge)


        return ret

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
