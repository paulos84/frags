$(document).ready(function(){
    let data = $('#proteins_button').data('proteins');
    $.each(data, function(key, val) {
        if (key > 10) {return false}
        let tds = '<td><button title=\"' + val.citation + ' PDB ID: ' + val.pdb_number + '\"><span class=\"glyphicon glyphicon-leaf\"></span></button></td><td><button type=\"button\" class=\"btn btn-primary enz-button\" value=\"' + val.pdb_number + '\">' + val.pdb_number + '</button></td><td style=\"padding-left:10px;\">' + val.notes + '</td>';
        $('#enz-table').append('<tr style="height:60px">' + tds + '</tr>');
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
    $('.dropdown-submenu a.action_parent').hover( function(event){
      event.preventDefault();
      // Avoid having the menu to close when clicking
      event.stopPropagation();
      // If a menu is already open we close it
      $('.dropdown-submenu').not('.mech_parent').removeClass('open');
      // opening the one you clicked on
      $(this).parent().addClass('open');
    });
    $('.mech_parent').hover( function(event){
      event.preventDefault();
      event.stopPropagation();
      $('.dropdown-submenu a.mech_parent').parent().removeClass('open');
      $(this).parent().addClass('open');
    });
    $('#chem_prop_button').click(function() {
        $('#chem_properties_table').toggle();
        $('#substructure_matches').hide();
        $('#user_data_table').hide();
        $('#3dmol-viewer').hide();
    })
    $('#user_data_button').click(function() {
        $('#user_data_table').toggle();
        $('#substructure_matches').hide();
        $('#chem_properties_table').hide();
        $('#3dmol-viewer').hide();
    })
    $('#substruct_match_button').click(function() {
        $('#substructure_matches').toggle();
        $('#chem_properties_table').hide();
        $('#user_data_table').hide();
        $('#3dmol-viewer').hide();
    })
    $('#three_d_button').click(function() {
        $('#3dmol-viewer').toggle();
        $('#chem_properties_table').hide();
        $('#substructure_matches').hide();
        $('#user_data_table').hide();
    })
});