Installation
============
The following is a brief guide to installing `BcForms`. The `Dockerfile <https://github.com/KarrLab/bpforms/blob/master/Dockerfile>`_ in the `BpForms` repository contains detailed instructions for how to install `BcForms` in Ubuntu Linux.

Prerequisites
--------------------------

First, install the third-party packages listed below.

* `ChemAxon Marvin <https://chemaxon.com/products/marvin>`_: optional to calculate major protonation and tautomerization states and draw molecules

    * `Java <https://www.java.com>`_ >= 1.8

* `Open Babel <http://openbabel.org>`_
* `Pip <https://pip.pypa.io>`_ >= 18.0
* `Python <https://www.python.org>`_ >= 3.6

To use ChemAxon Marvin to calculate major protonation and tautomerization states, set ``JAVA_HOME`` to the path to your Java virtual machine (JVM) and add Marvin to the Java class path::

   export JAVA_HOME=/usr/lib/jvm/default-java
   export CLASSPATH=$CLASSPATH:/opt/chemaxon/marvinsuite/lib/MarvinBeans.jar

Latest release From PyPI
---------------------------
Run the following command to install the latest release from PyPI.::

    pip install bcforms[all]

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub.::

    pip install git+https://github.com/KarrLab/pkg_utils.git#egg=pkg_utils
    pip install git+https://github.com/KarrLab/wc_utils.git#egg=wc_utils[chem, protonation]
    pip install git+https://github.com/KarrLab/bpforms.git#egg=bpforms
    pip install git+https://github.com/KarrLab/bcforms.git#egg=bcforms[all]

Installing the optional features
--------------------------------
To draw molecules, `BcForms` must be installed with the `[draw]` option:::

    pip install bcforms[draw]
    pip install git+https://github.com/KarrLab/bcforms.git#egg=bcforms[draw]


To install the REST API, `BcForms` must be installed with the `[rest_api option]`:::

    pip install bcforms[rest_api]
    pip install git+https://github.com/KarrLab/bcforms.git#egg=bcforms[rest_api]
