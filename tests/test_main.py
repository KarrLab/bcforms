""" Tests of bcforms command line interface (bcforms.__main__)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-06-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

from bcforms import __main__
import bcforms
import capturer
import mock
import unittest


class CliTestCase(unittest.TestCase):

    def test_cli(self):
        with mock.patch('sys.argv', ['bcforms', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: bcforms')

    def test_help(self):
        with self.assertRaises(SystemExit):
            with __main__.App(argv=['--help']) as app:
                app.run()

    def test_version(self):
        with __main__.App(argv=['-v']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), bcforms.__version__)
                self.assertEqual(captured.stderr.get_text(), '')

        with __main__.App(argv=['--version']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), bcforms.__version__)
                self.assertEqual(captured.stderr.get_text(), '')

    def test_validate(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['validate', '2 * abc_a + 3 * abc_b']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'Form is valid')
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['validate', 'bmp2_a + bmp2_a | crosslink: [l-bond-atom: bmp2_a(1)-362S1 | l-displaced-atom: bmp2_a(1)-362H1 | r-bond-atom: bmp2_a(2)-362S1 | r-displaced-atom: bmp2_a(2)-362H1]']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'Form is valid')
                self.assertEqual(captured.stderr.get_text(), '')

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['validate', 'HELLO']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['validate', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_c(1)-2O1 | l-displaced-atom: abc_d(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]']) as app:
                app.run()

    def test_get_formula(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-formula', 'abc_a + abc_b', '{abc_a:C5H10O, abc_b:C3H5O}']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'C8H15O2')
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-formula', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1+1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]', '{abc_a:C5H10O, abc_b:C3H5O}']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'C8H13O')
                self.assertEqual(captured.stderr.get_text(), '')

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['get-formula', 'HELLO', '{HELLO:C5H10O}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['get-formula', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_c(1)-2O1 | l-displaced-atom: abc_d(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]', '{abc_a:C5H10O, abc_b:C3H5O}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Cannot parse subunit_formulas:'):
            with __main__.App(argv=['get-formula', 'abc_a + abc_b', '{abc_a:C5H10O, abc_b}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Unable to calculate BcForm formula:'):
            with __main__.App(argv=['get-formula', 'abc_a + abc_b', '{abc_a:C5H10O}']) as app:
                app.run()

    def test_get_mol_wt(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-molwt', 'abc_a + abc_b', '{abc_a:86, abc_b:57}']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertAlmostEqual(float(captured.stdout.get_text()), 143, places=0)
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-molwt', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1+1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]', '{abc_a:86, abc_b:57}']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertAlmostEqual(float(captured.stdout.get_text()), 125, places=0)
                self.assertEqual(captured.stderr.get_text(), '')

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['get-molwt', 'HELLO', '{HELLO:86}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['get-molwt', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_c(1)-2O1 | l-displaced-atom: abc_d(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]', '{abc_a:86, abc_b:57}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Cannot parse subunit_mol_wts:'):
            with __main__.App(argv=['get-molwt', 'abc_a + abc_b', '{abc_a:86, abc_b}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Unable to calculate BcForm molecular weights:'):
            with __main__.App(argv=['get-molwt', 'abc_a + abc_b', '{abc_a:86}']) as app:
                app.run()

    def test_get_charge(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-charge', 'abc_a + abc_b', '{abc_a:-1, abc_b:0}']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(int(captured.stdout.get_text()), -1)
                self.assertEqual(captured.stderr.get_text(), '')

        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['get-charge', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_a(1)-2O1 | l-displaced-atom: abc_a(1)-2H1+1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]', '{abc_a:-1, abc_b:0}']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(int(captured.stdout.get_text()), -2)
                self.assertEqual(captured.stderr.get_text(), '')

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['get-charge', 'HELLO', '{HELLO:0}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
            with __main__.App(argv=['get-charge', 'abc_a + abc_b | crosslink: [l-bond-atom: abc_c(1)-2O1 | l-displaced-atom: abc_d(1)-2H1 | r-bond-atom: abc_b(1)-3C1 | r-displaced-atom: abc_b(1)-3H1 | r-displaced-atom: abc_b(1)-3O1]', '{abc_a:-1, abc_b:0}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Cannot parse subunit_charges:'):
            with __main__.App(argv=['get-charge', 'abc_a + abc_b', '{abc_a:-1, abc_b}']) as app:
                app.run()

        with self.assertRaisesRegex(SystemExit, '^Unable to calculate BcForm charges:'):
            with __main__.App(argv=['get-charge', 'abc_a + abc_b', '{abc_a:-1}']) as app:
                app.run()
