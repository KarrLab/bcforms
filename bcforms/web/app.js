$(document).foundation()

var num_rows=0
$('#add_subunit').click(function(){
   $('#table_subunits_dynamic').append(' \
        <tr id="row'+num_rows+'"> \
          <td> \
            Subunit Name: \
            <input type="text" id="subunit_name'+num_rows+'" name="subunit_name" placeholder="abc_a"/> \
            Structure encoding: \
            <input type="text" id="encoding'+num_rows+'" name="encoding" placeholder="bpforms.ProteinForm"/> \
            Structure string: \
            <input type="text" id="structure'+num_rows+'" name="structure" placeholder="AA"/> \
            Empirical Formula: \
            <input type="text" id="formula'+num_rows+'" name="formula"/> \
            Molecular weight: \
            <input type="text" id="mol_wt'+num_rows+'" name="mol_wt"/> \
            Charge: \
            <input type="text" id="charge'+num_rows+'" name="charge"/> \
          </td> \
          <td><button type="button" name="remove" id="'+num_rows+'" class="remove_subunit">X</button></td> \
        </tr>')
    num_rows++
})
$(document).on('click', '.remove_subunit', function(){
     var button_id = $(this).attr("id")
     $('#row'+button_id+'').remove()
})

$('#submit').click(function (evt) {
    bc_form = $('#bc_form_in').val().trim()

    if (bc_form == null || bc_form == '') {
        return
    }

    subunits = [];
    for (var i=0; i<num_rows; i++) {
        // subunit_name is required
        subunit_name = $('#subunit_name'+i+'').val().trim()
        if (subunit_name == null || subunit_name == '') {
            return
        }
        subunit = {'subunit_name':subunit_name}

        // other fields are optional
        encoding = $('#encoding'+i+'').val().trim()
        if (encoding != null && encoding != '') {
            subunit['encoding'] = encoding
        }
        structure = $('#structure'+i+'').val().trim()
        if (structure != null && structure != '') {
            subunit['structure'] = structure
        }
        formula = $('#formula'+i+'').val().trim()
        if (formula != null && formula != '') {
            subunit['formula'] = formula
        }
        mol_wt = $('#mol_wt'+i+'').val().trim()
        if (mol_wt != null && mol_wt != '') {
            subunit['mol_wt'] = parseFloat(mol_wt)
        }
        charge = $('#charge'+i+'').val().trim()
        if (charge != null && charge != '') {
            subunit['charge'] = parseInt(charge)
        }
        subunits.push(subunit)
    }

    data = {
        'form': bc_form,
        'subunits': subunits
    };

    console.log(JSON.stringify(data))

    // $.ajax({
    //   type: 'post',
    //   url: '/api/bcform/',
    //   data: JSON.stringify(data),
    //   contentType : 'application/json',
    //   dataType: 'json',
    //   success: set_properties
    // })

})
// 
// set_properties = function(data, status, jqXHR) {
//     $("#output_test").val(status)
// }
