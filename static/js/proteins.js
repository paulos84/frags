$(document).ready(function(){
    const data = $('#proteins_button').data('proteins');
    $.each(data, function(key, val) {
        let tds = '<td><button title=\"' + val.citation + ' PDB ID: ' + val.pdb_number + '\"><span class=\"glyphicon glyphicon-leaf\"></span></button></td><td><button type=\"button\" class=\"btn btn-primary enz-button\" value=\"' + val.pdb_number + '\">' + val.pdb_number + '</button></td><td>' + val.notes + '</td>';
        $('#enz-table').append('<tr style="height:60px">' + tds + '</tr>');
    });
    $('body').on('click', 'button.enz-button', function(){
        $("#enz-loading").show();
        $("html, body").scrollTop($("#enz-loading").offset().top - 64);
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
});