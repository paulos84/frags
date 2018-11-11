$(document).ready(function(){

    var height = $('#cores-column').height()
    var height2 = $('#form_sidebar').height()
    $('.sidebar-nav').height(height + height2)

    $("#id_property_choice").change(function () {
      var choice = $("#id_property_choice").val();
      $(".chem_data_plot").hide();
      $("#" + choice).show();
      $("#sd_" + choice).show();
    });
});