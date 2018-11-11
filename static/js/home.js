$(document).ready(function() {

  $('.dropdown-submenu a.action_parent').hover( function(event){
      event.preventDefault();
      // Avoid having the menu to close when clicking
      event.stopPropagation();
      // If a menu is already open we close it
      $('.dropdown-submenu').removeClass('open');
      // opening the one you clicked on
      $(this).parent().addClass('open');
  });

});