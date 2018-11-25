$(document).ready(function(){
    var height = $('#cores-column').height()
    var height2 = $('#form_sidebar').height()
    $('.col-sm-2').height(height + height2)
    $('.dropdown-submenu a.action_parent').hover( function(event){
      event.preventDefault();
      event.stopPropagation();
      $('.dropdown-submenu').not('.mech_parent').removeClass('open');
      $(this).parent().addClass('open');
    });
    $('.mech_parent').hover( function(event){
      event.preventDefault();
      event.stopPropagation();
      $('.dropdown-submenu a.mech_parent').parent().removeClass('open');
      $(this).parent().addClass('open');
    });
    $("#id_property_choice").change(function () {
      var choice = $("#id_property_choice").val();
      $(".chem_data_plot").hide();
      $("#" + choice).show();
      $("#sd_" + choice).show();
    });
});