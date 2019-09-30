# Import libraries
import bcforms
import bpforms

# Create complexes from their string representations
form_1 = bcforms.BcForm().from_str('2 * subunit_a + 3 * subunit_b')
form_1.set_subunit_attribute('subunit_a', 'structure',
    bpforms.ProteinForm().from_str('CAAAAAAAA'))
form_1.set_subunit_attribute('subunit_b', 'structure',
    bpforms.ProteinForm().from_str('AAAAAAAAC'))

form_2 = bcforms.BcForm().from_str(
    '2 * subunit_a'
    '| x-link: [type: disulfide | l: subunit_a(1)-1 | r: subunit_a(2)-1]')
form_2.set_subunit_attribute('subunit_a', 'structure',
    bpforms.ProteinForm().from_str('CAAAAAAAA'))

# Create complexes programmatically
form_1_b = bcforms.BcForm()
form_1_b.subunits.append(bcforms.core.Subunit('subunit_a', 2, 
    bpforms.ProteinForm().from_str('CAAAAAAAA')))
form_1_b.subunits.append(bcforms.core.Subunit('subunit_b', 3, 
    bpforms.ProteinForm().from_str('AAAAAAAAC')))

form_2_b = bcforms.BcForm()
subunit = bcforms.core.Subunit('subunit_a', 2, 
    bpforms.ProteinForm().from_str('CAAAAAAAA'))
form_2_b.subunits.append(subunit)
form_2_b.crosslinks.append(bcforms.core.OntologyCrosslink(
    'disulfide', 'subunit_a', 1, 'subunit_a', 1, 1, 2))

# Get properties of polymers
form_1.subunits[0].id # >> subunit_a
form_2.crosslinks[0] # >> <bcforms.core.OntologyCrosslink at 0x7f89c56949e8>

# Get the string representation of a complex
str(form_1_b) # >> 2 * subunit_a + 3 * subunit_b

# Check equality of polymers
form_1_b.is_equal(form_1) # >> True

# Calculate properties of a polymer
form_1.get_structure()[0] # >> <openbabel.OBMol>
form_1.export('smiles') # >> C(=O)([C@@H]([NH3+])CS)N...
str(form_1.get_formula()) # >> C135H240N45O50S5
form_1.get_charge() # >> 5
