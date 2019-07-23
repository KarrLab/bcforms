$(document).foundation()

var num_rows=0
$('#add_subunit').click(function(){
   $('#table_subunits_dynamic').append(' \
        <tr id="row'+num_rows+'"> \
          <td> \
            Subunit Name: \
            <input type="text" id="subunit_name'+num_rows+'" name="subunit_name" placeholder="abc_a"/> \
            <label>Select type of known information \
              <select class= "subunit_info_select" id="subunit_info_type_'+num_rows+'"> \
                <option value=1>Physical properties</option> \
                <option value=0 selected>Structure</option> \
              </select> \
            </label> \
            <div id="subunit_info'+num_rows+'"> \
              Structure encoding: \
              <input type="text" id="encoding'+num_rows+'" name="encoding" placeholder="bpforms.ProteinForm"/> \
              Structure string: \
              <input type="text" id="structure'+num_rows+'" name="structure" placeholder="AA"/> \
            </div> \
          </td> \
          <td><button type="button" name="remove" id="'+num_rows+'" class="remove_subunit">X</button></td> \
        </tr>')
    num_rows++
})

$(document).on('click', '.remove_subunit', function(){
     var button_id = $(this).attr("id")
     $('#row'+button_id+'').remove()
})

$(document).on('change', '.subunit_info_select', function(){
    var select_id = $(this).attr("id").substring(18)
    if ($('#subunit_info_type_'+select_id).val() == 0) {
        $('#subunit_info'+select_id).html(' \
            Structure encoding: \
            <input type="text" id="encoding'+select_id+'" name="encoding" placeholder="bpforms.ProteinForm"/> \
            Structure string: \
            <input type="text" id="structure'+select_id+'" name="structure" placeholder="AA"/>')
    }
    else {
        $('#subunit_info'+select_id).html(' \
            Empirical Formula: \
            <input type="text" id="formula'+select_id+'" name="formula"/> \
            Molecular weight: \
            <input type="text" id="mol_wt'+select_id+'" name="mol_wt"/> \
            Charge: \
            <input type="text" id="charge'+select_id+'" name="charge"/>')
    }
})


$('#submit').click(function (evt) {
    bc_form = $('#bc_form_in').val().trim()

    if (bc_form == null || bc_form == '') {
        return
    }

    subunits = [];
    for (var i=0; i<num_rows; i++) {
        if ($('#row'+i).length) {
            // subunit_name is required
            subunit_name = $('#subunit_name'+i+'').val().trim()
            if (subunit_name == null || subunit_name == '') {
                return
            }
            subunit = {'subunit_name':subunit_name}

            // other fields are optional
            if ($('#encoding'+i+'').length && typeof $('#encoding'+i+'').val() !== 'undefined') {
                encoding = $('#encoding'+i+'').val().trim()
                if (encoding != null && encoding != '') {
                    subunit['encoding'] = encoding
                }
            }
            if ($('#structure'+i+'').length && typeof $('#structure'+i+'').val() !== 'undefined' ) {
                structure = $('#structure'+i+'').val().trim()
                if (structure != null && structure != '') {
                    subunit['structure'] = structure
                }
            }
            if ($('#formula'+i+'').length && typeof $('#formula'+i+'').val() !== 'undefined') {
                formula = $('#formula'+i+'').val().trim()
                if (formula != null && formula != '') {
                    subunit['formula'] = formula
                }
            }
            if ($('#mol_wt'+i+'').length && typeof $('#mol_wt'+i+'').val() !== 'undefined') {
                mol_wt = $('#mol_wt'+i+'').val().trim()
                if (mol_wt != null && mol_wt != '') {
                    subunit['mol_wt'] = parseFloat(mol_wt)
                }
            }
            if ($('#charge'+i+'').length && typeof $('#charge'+i+'').val() !== 'undefined') {
                charge = $('#charge'+i+'').val().trim()
                if (charge != null && charge != '') {
                    subunit['charge'] = parseInt(charge)
                }
            }

            subunits.push(subunit)
        }
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
