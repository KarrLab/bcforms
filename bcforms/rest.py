""" REST JSON API

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-07-03
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bcforms
import bcforms.core
import bpforms
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

# define model

# if encoding, structure defined -> ignore formula, mol_wt, charge, and define them based on structure
# if neither encoding, structure set and formula is defined -> ignore mol_wt, and define mol_wt based on formula
subunit_fields = {}
subunit_fields['subunit_name'] = flask_restplus.fields.String(required=True, title='Subunit name', example='abc_a')
# encoding can be smiles, bpforms.ProteinForm, bpforms.DnaForm, bpforms.RnaForm
subunit_fields['encoding'] = flask_restplus.fields.String(required=False, title='Structure encoding', example='bpforms.ProteinForm')
subunit_fields['structure'] = flask_restplus.fields.String(required=False, title='Structure string', example='AAA')
subunit_fields['formula'] = flask_restplus.fields.String(required=False, title='Empirical formula', example='C5H10O')
subunit_fields['mol_wt'] = flask_restplus.fields.Float(required=False, title='Molecular weight', example=86.0)
subunit_fields['charge'] = flask_restplus.fields.Integer(required=False, title='Total charge', example=0)

bcform_fields = {}
bcform_fields['form'] = flask_restplus.fields.String(required=True, title='BcForm', description='input biocomplex form', example='2 * abc_a + 3 * abc_b')
bcform_fields['subunits'] = flask_restplus.fields.List(flask_restplus.fields.Nested(bcform_ns.model('Subunit',subunit_fields)), example=[
    {
      "subunit_name": "abc_a",
      "encoding": "bpforms.ProteinForm",
      "structure": "AAA"
    },
    {
      "subunit_name": "abc_b",
      "encoding": "bpforms.ProteinForm",
      "structure": "MM"
    }
  ])

bcforms_model = bcform_ns.model('BcForm', bcform_fields)


@bcform_ns.route("/")
class Bcform(flask_restplus.Resource):

    @bcform_ns.expect(bcforms_model, validate=True)
    def post(self):
        ret = {}

        args = bcform_ns.payload

        # print(args)

        # get arguments
        form = args['form']
        arg_subunits = args.get('subunits', None)

        # validate form
        try:
            bc_form = bcforms.core.BcForm().from_str(form)
        except Exception as error:
            flask_restplus.abort(400, 'Form is invalid', errors={'form': str(error)})

        errors = bc_form.validate()
        if errors:
            flask_restplus.abort(400, 'Form is invalid', errors={'form': '. '.join(errors)})

        # validate input subunit properties
        if arg_subunits is not None:
            for subunit in arg_subunits:

                # check if name is in the form
                subunit_id = subunit['subunit_name']
                if subunit_id in [subunit.id for subunit in bc_form.subunits]:

                    # check if encoding and structure are present at the same time
                    if ('encoding' in subunit) and ('structure' in subunit):
                        # if encoding and structure both present, check if encoding is known
                        encoding = subunit['encoding'].strip()
                        if encoding == 'bpforms.ProteinForm':
                            try:
                                subunit_structure = bpforms.ProteinForm().from_str(subunit['structure'])
                                bc_form.set_subunit_attribute(subunit_id, 'structure', subunit_structure)
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse bpforms.ProteinForm', errors={'structure': str(error)})
                        elif encoding == 'bpforms.DnaForm':
                            try:
                                subunit_structure = bpforms.DnaForm().from_str(subunit['structure'])
                                bc_form.set_subunit_attribute(subunit_id, 'structure', subunit_structure)
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse bpforms.DnaForm', errors={'structure': str(error)})
                        elif encoding == 'bpforms.RnaForm':
                            try:
                                subunit_structure = bpforms.RnaForm().from_str(subunit['structure'])
                                bc_form.set_subunit_attribute(subunit_id, 'structure', subunit_structure)
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse bpforms.RnaForm', errors={'structure': str(error)})
                        elif encoding == 'smiles' or encoding == 'SMILES' or encoding == 'smi' or encoding == 'SMI':
                            try:
                                bc_form.set_subunit_attribute(subunit_id, 'structure', subunit['structure'])
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse SMILES string', errors={'structure': str(error)})

                    # else if one is present but not the other, report error
                    elif ('encoding' in subunit) ^ ('structure' in subunit):
                        flask_restplus.abort(400, 'One of encoding and structure is present but not both')

                    # when neither encoding nor structure is present
                    else:
                        # check formula
                        if 'formula' in subunit:
                            try:
                                bc_form.set_subunit_attribute(subunit_id, 'formula', subunit['formula'])
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse formula', errors={'formula': str(error)})
                        elif 'mol_wt' in subunit:
                            try:
                                bc_form.set_subunit_attribute(subunit_id, 'mol_wt', subunit['mol_wt'])
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse mol_wt', errors={'mol_wt': str(error)})

                        # check charge
                        if 'charge' in subunit:
                            try:
                                bc_form.set_subunit_attribute(subunit_id, 'charge', subunit['charge'])
                            except Exception as error:
                                flask_restplus.abort(400, 'Unable to parse charge', errors={'charge': str(error)})

                else:
                    flask_restplus.abort(400, 'Subunit name not in BcForm')


        ret['form'] = str(bc_form)

        try:
            ret['structure'] = bc_form.export()
        except Exception:
            pass

        try:
            ret['formula'] = str(bc_form.get_formula())
        except Exception:
            pass

        try:
            ret['mol_wt'] = bc_form.get_mol_wt()
        except Exception:
            pass

        try:
            ret['charge'] = bc_form.get_charge()
        except Exception:
            pass

        return ret

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
