""" Test examples of BcForms

:Author: Mike Zheng <xzheng20@colby.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-07-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

import glob
import itertools
import json
import nbconvert.preprocessors
import nbformat
import os
import shutil
import sys
import tempfile
import unittest

@unittest.skipIf(os.getenv('CIRCLECI', '0') in ['1', 'true'], 'Jupyter server not setup in CircleCI')
class ExamplesTestCase(unittest.TestCase):
    TIMEOUT = 600

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, 'examples')

    @classmethod
    def tearDownClass(cls):
        sys.path.remove('examples')

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_jupyter(self):
        for filename in itertools.chain(glob.glob('examples/*.ipynb'), glob.glob('examples/**/*.ipynb')):
            with open(filename) as file:
                version = json.load(file)['nbformat']
            with open(filename) as file:
                notebook = nbformat.read(file, as_version=version)
            execute_preprocessor = nbconvert.preprocessors.ExecutePreprocessor(timeout=self.TIMEOUT)
            execute_preprocessor.preprocess(notebook, {'metadata': {'path': 'examples/'}})
