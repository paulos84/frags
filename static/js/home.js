$(document).ready(function() {
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
  $('#development_table').dataTable({
        paging: true,
        bAutoWidth: false ,
        aoColumns : [
            { sWidth: '5%' }, 
			{ sWidth: '5%' }, 
			{ sWidth: '12%' }, 
			{ sWidth: '18%' },
			{ sWidth: '12%' },
            { sWidth: '6%' },
			{ sWidth: '42%' }
        ],
        dom: 'Bfrtip',
        buttons: [
            {
			 text: 'Show current trials',
			 action: function () {
				$('.study-col').toggle();
				$('.dt-buttons').children().first().text(function(i, text){
					  return text === "Hide current trials" ? "Show current trials" : "Hide current trials";
				  });
				}
			},
            {
			 extend: 'collection',
             text: 'Select company',
             buttons: companies
			 },
        ]
    });
});
