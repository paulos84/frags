
var checkedCheckboxes = [];

$(document).ready(function(){

    var height = $('#compounds_list').height();
    $('.col-sm-2').height(height + 800);

    var table_height = $('#enz-table').height();
    var viewer_height = $('#mol-viewer').height();
    if (table_height > viewer_height) {
        $('#mol-viewer').height(table_height);
    } else if (table_height < viewer_height) {
        $('#enz-table').height(viewer_height);
    }

    $("#show_enz").click(function () {
      $(".enz_view").show();
      $('.col-sm-2').height(height + mol_height + plot_height + 550);
      $(this).hide();
    });

    $("#chem_compare").click(function () {
      $(".compound_checkbox").show();
      $("#submit_selection").show();
      $("#chem_compare").hide();
    });

    $("#id_property_choice").change(function () {
        var choice = $("#id_property_choice").val();
        $(".chem_data_plot").hide();
        $("#" + choice).show();
    });

    $("#show_activities").click(function () {
        $(this).text(function(i, text){
        return text === "Hide mechanisms" ? "Show mechanisms" : "Hide mechanisms";
      })
    });

    $(".enz-button").click(function() {
        $("#enz-loading").show();
        let pdb_no = $(this).val();
        let element = $('#container-01');
        let config = {};
        let viewer = $3Dmol.createViewer( element, config );
        $3Dmol.download("pdb:" + pdb_no,viewer,{},function(){
              viewer.setStyle({chain: 'A', within:{distance: 10, sel:{chain: 'B'}}}, {sphere:{}});
              viewer.setStyle({chain:'B'},{cartoon:{color:'spectrum'}});
              viewer.setStyle({chain:'B',invert:true},{cartoon:{}});
              viewer.render();
              $("#enz-loading").hide();
            });
    });

    $('#submit_button').click(function(){
        checkedCheckboxes = []
        $("input").each(function(){
            console.log('????')
            if ($(this).is(':checked')) {
                checkedCheckboxes.push($(this).val())
            }
        });
    });
});

