Command line interface
----------------------

The command line interface provides five functions to easily manipulate `BcForms`-encoded descriptions of complexes.

* **Get help with the `BcForms` command line interface.** The following commands return inline help information about the command line interface::

    bcforms
    bcforms -h
    bcforms --help

* **Validate a `BcForms`-encoded description of a form of a complex.** The following command can be used to verify if description of a complex is syntactically and semantically valid. The command line interface will print any errors to the standard error::

    bcforms validate <bcform>

    bcforms validate '2 * a + 3 * b'
    # Form is valid

* **Calculate the formula of a complex.** The following command can be used to calculate the formula of a complex::

    bcforms get-formula --help

    bcforms get-formula <bcform> <dictionary of formulae of subunits>

    # Calculate the formula of a complex
    bcforms get-formula '2 * a + 3 * b' '{a: CHO, b: C2H2O2}'
    # C8H8O8

* **Calculate the charge of a complex.** The following command can be used to calculate the charge of a complex::

    bcforms get-charge --help

    bcforms get-charge <bcform> <dictionary of charges of subunits>

    # Calculate the charge of a complex
    bcforms get-charge '2 * a + 3 * b' '{a: 1, b: 2}'
    # 8

* **Calculate the molecular weight of a complex.** The following command can be used to calculate the molecular weight of a complex::

    bcforms get-molwt --help
    bcforms get-molwt <bcform> <dictionary of molecular weights of subunits>

    # Calculate the molecular weight of a complex
    bcforms get-molwt '2 * a + 3 * b' '{a: 1, b: 2}'
    # 8
