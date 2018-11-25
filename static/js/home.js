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
});