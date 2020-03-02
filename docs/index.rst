`BcForms` documentation
=======================

`BcForms` is a toolkit for concretely describing the molecular structure (atoms and bonds) of macromolecular complexes, including non-canonical monomeric forms, circular topologies, and crosslinks. `BcForms` was developed to help describe the semantic meaning of whole-cell computational models (see `https://wholecell.org <https://www.wholecell.org>`_).

`BcForms` includes a grammar for describing forms of macromolecular complexes composed of DNA, RNA, protein, and small molecular subunits and crosslinks between the subunits. The DNA, RNA, and protein subunits can be described using `BpForms <https://www.bpforms.org>`_ and the small molecule subunits can be described using SMILES. `BcForms` also includes four software tools for verifying descriptions of complexes and calculating physical properties of complexes such as their molecular structure, formula, molecular weight, and charge: this website, a `JSON REST API <https://www.bcforms.org/api/>`_, a command line interface, and a Python API. `BcForms` is available open-source under the MIT license.

Contents
--------

.. toctree::
   :maxdepth: 3
   :numbered:

   installation.rst
   cli.rst
   API documentation <source/modules.rst>
   about.rst
   references.rst