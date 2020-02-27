[//]: # ( [![PyPI package](https://img.shields.io/pypi/v/bcforms.svg)](https://pypi.python.org/pypi/bcforms) )
[![Documentation](https://readthedocs.org/projects/bcforms/badge/?version=latest)](https://docs.karrlab.org/bcforms)
[![Test results](https://circleci.com/gh/KarrLab/bcforms.svg?style=shield)](https://circleci.com/gh/KarrLab/bcforms)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/bcforms/badge.svg)](https://coveralls.io/github/KarrLab/bcforms)
[![Code analysis](https://api.codeclimate.com/v1/badges/c8f15ac1e50c27ca44cc/maintainability)](https://codeclimate.com/github/KarrLab/bcforms)
[![License](https://img.shields.io/github/license/KarrLab/bcforms.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/bcforms/README.md?pixel)

# `BcForms`: toolkit for concretely describing macromolecular complexes

`BcForms` is a toolkit for concretely describing the molecular structure (atoms and bonds) of macromolecular complexes, including non-canonical monomeric forms, circular topologies, and crosslinks.

`BcForms` encompasses five tools:
* A grammar for concretely describing subunit composition and interchain crosslinks of biocomplexes. See the [documentation](https://docs.karrlab.org/bcforms/) for more information. For example, the following text represents a homodimer complex with a disulfide bond between the first monomer in the first subunit and the first monomer in the second subunit.
    ```
    complex: 2 * sub_c | x-link: [
        l-bond-atom: sub_c(1)-1S11 |
        l-displaced-atom: sub_c(1)-1H11 |
        r-bond-atom: sub_c(2)-1S11 |
        r-displaced-atom: sub_c(2)-1H11
      ]
    ```

    This concrete representation enables the `BcForms` software tools to calculate properties of biocomplexes when the subunits are concretely defined.

* Tools for calculating properties of biocomplexes including their structure, chemical formulae, molecular weights, and charges.
  * A web app: [https://bcforms.org](https://bcforms.org)
  * A JSON REST API: [https://bcforms.org/api](https://bcforms.org/api)
  * A command line interface. See the [documentation](https://docs.karrlab.org/bcforms/master/0.0.1/) for more information.
  * A Python API. See the [documentation](https://docs.karrlab.org/bcforms/master/0.0.1/) for more information.

## Installation
The following is a brief guide to installing `BcForms`. The [Dockerfile](https://github.com/KarrLab/bpforms/blob/master/Dockerfile) in the `BpForms` repository contains detailed instructions for how to install `BcForms` in Ubuntu Linux.

1. Install dependencies

   * [ChemAxon Marvin](https://chemaxon.com/products/marvin): optional to calculate major protonation and tautomerization states and draw molecules
     * [Java](https://www.java.com) >= 1.8
   * [Open Babel](http://openbabel.org)
   * [Pip](https://pip.pypa.io) >= 19.0
   * [Python](https://www.python.org) >= 3.6

2. Install this package

    * Install the latest release from PyPI
      ```
      pip install bcforms[all]
      ```

    * Install the latest revision from GitHub
      ```
      pip install git+https://github.com/KarrLab/pkg_utils.git#egg=pkg_utils
      pip install git+https://github.com/KarrLab/wc_utils.git#egg=wc_utils[chem, protonation]
      pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms
      pip install git+https://github.com/KarrLab/bcforms.git#egg=bcforms[all]
      ```

    * To draw molecules, `BcForms` must be installed with the `[draw]` option:
      ```
      pip install bcforms[draw]
      pip install git+https://github.com/KarrLab/bcforms.git#egg=bcforms[draw]
      ```

    * To install the REST API, `BcForms` must be installed with the `[rest_api option]`:
      ```
      pip install bcforms[rest_api]
      pip install git+https://github.com/KarrLab/bcforms.git#egg=bcforms[rest_api]
      ```

## Documentation
Please see the [API documentation](https://docs.karrlab.org/bcforms).

## License
The package is released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

* [Jonathan Karr](https://www.karrlab.org)
* [Xiaoyue Zheng](https://www.linkedin.com/in/xiaoyue-zheng/)

## Questions and comments
Please contact the [Karr Lab](mailto:info@karrlab.org) with any questions or comments.
