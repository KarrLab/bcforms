$(document).foundation()

$('#submit').click(function (evt) {
    bc_form = $("#bc_form_in").val().trim()

    if (bc_form == null || bc_form == '') {
        return;
    }

    data = {
        'form': bc_form
    }

    $.ajax({
      type: 'post',
      url: '/api/bcform/',
      data: JSON.stringify(data),
      contentType : 'application/json',
      dataType: 'json',
      success: set_properties
    })

})

set_properties = function(data, status, jqXHR) {
    $("#output_test").val(status)
}
