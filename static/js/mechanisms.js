var checkedCheckboxes = [];
$(document).ready(function(){
    var height = $('#mechanisms-column').height()
    $('.col-sm-2').height(height);
        var height = $('#cores-column').height()
    var height2 = $('#form_sidebar').height()
    setProteinHeight()
    $("#show-enz").click(function() {
      $(this).text(function(i, text){
          return text === "Compound entry counts" ? "View target protein structures" : "Compound entry counts";
      });
      $("#error_text").hide();
      if ($('.enz-panel-button-container').is(':visible')) {
            $('.mech-count').show();
            $('.enz-panel-button-container').hide();
      } else {
       $('.mech-count').hide();
       $('.enz-panel-button-container').show();
      }
      $('.mech-lit-button-container').hide();
      $("#info_text").toggle();
    });
    $("#show_lit").click(function() {
      $('.mech-count').hide();
      $('.enz-panel-button-container').hide();
      $('.mech-lit-button-container').show();
    });
    $('body').on('click', 'button.enz-button', function(){
        $("#enz-loading").show();
        let pdb_no = $(this).val();
        let element = $('#container-01');
        let config = {};
        let viewer = $3Dmol.createViewer( element, config );
        $3Dmol.download("pdb:" + pdb_no,viewer,{},function(){
              viewer.setStyle({chain: 'B', within:{distance: 10, sel:{chain: 'A'}}}, {sphere:{}});
              viewer.setStyle({chain:'A'},{cartoon:{color:'spectrum'}});
              viewer.setStyle({chain:'A',invert:true},{cartoon:{}});
              viewer.render();
              $("#enz-loading").hide();
              setProteinHeight();
            });
    });
    $("#chem_prop").click(function(){
      $("#chem_prop").hide()
      $(".compare_prop").show();
      $("#submit_button").show();
      $("#sel_info_text").show();
    });
    $(".compare_prop").click(function(){
        $("#error_text").hide();
    });
    $("input:checkbox[name='selected_mechanisms']").click(function(){
        $("#sel_info_text").hide();
    })
    $('#submit_button').click(function(){
        checkedCheckboxes = []
        $("input:checkbox[name='selected_mechanisms']").each(function(){
            if ($(this).is(':checked')) {
                checkedCheckboxes.push($(this).val())
            }
        });
    });
    $("#id_property_choice").change(function (){
      var choice = $("#id_property_choice").val();
      $(".chem_data_plot").hide();
      $("#" + choice).show();
      $("#sd_" + choice).show();
    });
});
function formValidate(){
    if (checkedCheckboxes.length < 1 || checkedCheckboxes.length > 20) {
        $("#error_text").show();
        return false;
    };
};
function setProteinHeight() {
    var table_height = $('#enz-table').height();
    var viewer_height = $('#mol-viewer').height();
    if (table_height > viewer_height) {
        $('#mol-viewer').height(table_height);
    } else if (table_height < viewer_height) {
        $('#enz-table').height(viewer_height);
    }
}