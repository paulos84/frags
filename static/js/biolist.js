var checkedCheckboxes = [];

$(document).ready(function(){
    $('#proteins_button').click(function(){
        var data = $(this).data('proteins');
        $.each(data, function(key, val) {
            let tds = '<td><button title=\"' + val.citation + ' PDB ID: ' + val.pdb_number + '\"><span class=\"glyphicon glyphicon-leaf\"></span></button></td><td><button type=\"button\" class=\"btn btn-primary enz-button\" value=\"' + val.pdb_number + '\">' + val.pdb_number + '</button></td><td>' + val.notes + '</td>';
            $('#enz-table').append('<tr style="height:60px">' + tds + '</tr>');
        });
        $(".enz_view").show();
        $(this).hide();
        setProteinHeight();
        $("html, body").scrollTop(0);
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
    $("#chem_compare").click(function () {
      $(".compound_checkbox").show();
      $("#submit_selection").show();
      $("#chem_compare").hide();
    });
    $('#submit_button').click(function(){
        checkedCheckboxes = []
        $("input:checkbox[name='selected_bioactives']").each(function(){
            if ($(this).is(':checked')) {
                checkedCheckboxes.push($(this).val())
            }
        });
    });
});
function setProteinHeight() {
    var table_height = $('#enz-table').height();
    var viewer_height = $('#mol-viewer').height();
    if (table_height > viewer_height) {
        $('#mol-viewer').height(table_height);
    } else if (table_height < viewer_height) {
        $('#enz-table').height(viewer_height);
    }
    $('.col-sm-2').height($('#compounds_list').height() + viewer_height);
}
function formValidate(){
    if (checkedCheckboxes.length < 1 || checkedCheckboxes.length > 20) {
        $("#error_text").show();
        return false;
    };
};
