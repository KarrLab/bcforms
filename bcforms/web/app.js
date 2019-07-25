$(document).foundation()

var num_rows=0
$('#add_subunit').click(function(){
   $('#table_subunits_dynamic').append(' \
        <tr id="row'+num_rows+'"> \
          <td> \
            <label>Name:</label> \
            <input type="text" id="name'+num_rows+'" name="name" placeholder="abc_a"/> \
            \
            <label>Structural information</label> \
              <select class= "subunit_info_select" id="subunit_info_type_'+num_rows+'"> \
                <option value=1>Formula, molecular weight, and/or charge</option> \
                <option value=0 selected>Molecular structure</option> \
              </select> \
            \
            <div id="subunit_info'+num_rows+'"> \
              <label>Structure format</label> \
              <select name="encoding" class="encoding" id="encoding'+num_rows+'"> \
                <option value=3>SMILES</option> \
                <option value=2>bpforms.RnaForm</option> \
                <option value=1>bpforms.DnaForm</option> \
                <option value=0 selected>bpforms.ProteinForm</option> \
              </select> \
              \
              <label>Structure</label> \
              <input type="text" id="structure'+num_rows+'" name="structure" placeholder="AA"/> \
            </div> \
          </td> \
          <td><button type="button" name="remove" id="'+num_rows+'" class="alert button remove_subunit">Delete</button></td> \
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
            <label>Structure format</label> \
            <select name="encoding" class="encoding" id="encoding'+num_rows+'"> \
              <option value=3>SMILES</option> \
              <option value=2>bpforms.RnaForm</option> \
              <option value=1>bpforms.DnaForm</option> \
              <option value=0 selected>bpforms.ProteinForm</option> \
            </select> \
            \
            <label>Structure</label> \
            <input type="text" id="structure'+select_id+'" name="structure" placeholder="AA"/>')
    }
    else {
        $('#subunit_info'+select_id).html(' \
            <label>Formula</label> \
            <input type="text" class="formula_input" id="formula'+select_id+'" name="formula"/> \
            <label>Molecular weight</label> \
            <input type="text" id="mol_wt'+select_id+'" name="mol_wt"/> \
            <label>Charge</label> \
            <input type="text" id="charge'+select_id+'" name="charge"/>')
    }
})

$(document).on('change', '.formula_input', function(){
    var formula_id = $(this).attr("id").substring(7)
    if ($(this).val()) {
    	$('#mol_wt'+formula_id).val('')
        $('#mol_wt'+formula_id).prop("disabled", true)
    }
    else {
        $('#mol_wt'+formula_id).prop("disabled", false)
    }
})


$('#submit').click(function (evt) {
    bc_form = $('#bc_form_in').val().trim()

    if (bc_form == null || bc_form == '') {
        return
    }

    subunits = []
    for (var i=0; i<num_rows; i++) {
        if ($('#row'+i).length) {
            // name is required
            name = $('#name'+i+'').val().trim()
            if (name == null || name == '') {
                return
            }
            subunit = {'name': name}

            // other fields are optional
            if ($('#encoding'+i+'').length && typeof $('#encoding'+i+'').val() !== 'undefined') {
                encoding = $('#encoding'+i+'').val()
                if (encoding == 0) {
                    subunit['encoding'] = 'bpforms.ProteinForm'
                }
                else if (encoding == 1) {
                    subunit['encoding'] = 'bpforms.DnaForm'
                }
                else if (encoding == 2) {
                    subunit['encoding'] = 'bpforms.RnaForm'
                }
                else if (encoding == 3) {
                    subunit['encoding'] = 'SMILES'
                }
                else {
                    return
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
                mol_wt = parseFloat($('#mol_wt'+i+'').val().trim())
                if (mol_wt != null && (!isNaN(mol_wt))) {
                    subunit['mol_wt'] = mol_wt
                }
            }
            if ($('#charge'+i+'').length && typeof $('#charge'+i+'').val() !== 'undefined') {
                charge = parseInt($('#charge'+i+'').val().trim())
                if (charge != null && (!isNaN(charge))) {
                    subunit['charge'] = charge
                }
            }

            // validate subunit
            // if encoding is set, check if structure is set
            if (('encoding' in subunit) && (!('structure' in subunit))) {
                return
            }

            subunits.push(subunit)
        }
    }

    data = {
        'form': bc_form,
        'subunits': subunits
    }

//     console.log(JSON.stringify(data))

    $.ajax({
      type: 'post',
      url: '/api/bcform/',
      data: JSON.stringify(data),
      contentType : 'application/json',
      dataType: 'json',
      success: set_properties
    })
    .fail(display_error);

})

display_error = function( jqXHR, textStatus, errorThrown ) {
    error = '<b>' + jqXHR['responseJSON']['message'] + '</b>'
    if ('errors' in jqXHR['responseJSON']) {
        error += '<ul>'
        for (field in jqXHR['responseJSON']['errors'])
            error += '<li>'
            if (field != '')
                error += '<span style="text-decoration: underline;">' + field + '</span>: '
            error += jqXHR['responseJSON']['errors'][field]
            error += '</li>'
        error += '</ul>'
    }
    $("#errors").html(error)
    $("#errors").css('padding-bottom', '16px')
}

set_properties = function(data, status, jqXHR) {

    $("#errors").html('')
    $("#errors").css('padding-bottom', '0px')

	// clear everything
	$("#out_bcform").val('')
	$("#out_structure").val('')
	$("#out_formula").val('')
	$("#out_mol_wt").val('')
	$("#out_charge").val('')

	// write form
	form = data['form']
	$("#out_bcform").val(form)

	// write structure
	if ('structure' in data) {
        if (data['structure'] != null) {
    		structure = data['structure']
    		$("#out_structure").val(structure)
        }
	}

	// write formula
	if ('formula' in data) {
		formula = data['formula']
		$("#out_formula").val(formula)
	}

	// write mol_wt
	if ('mol_wt' in data) {
		mol_wt = data['mol_wt']
		$("#out_mol_wt").val(mol_wt)
	}

	// write charge
	if ('charge' in data) {
		charge = data['charge']
		$("#out_charge").val(charge)
	}

    if ('warnings' in data && data['warnings'] != null) {
        warnings = 'Warning: ' + data['warnings'] + ' Please download <i>BcForms</i> to calculate structure of longer biocomplexes.'
        $("#warnings").html(warnings)
        $("#warnings").css('padding-bottom', '16px')
    } else {
        warnings = ''
        $("#warnings").html(warnings)
        $("#warnings").css('padding-bottom', '0px')
    }


}
