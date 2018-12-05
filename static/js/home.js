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
            { sWidth: '5%' }, { sWidth: '5%' }, { sWidth: '12%' }, { sWidth: '18%' }, { sWidth: '12%' },
            { sWidth: '6%' }, { sWidth: '42%' }
        ],
        dom: 'Bfrtip',
        buttons: [
            {text: 'Show study info',
             action: function () {
                $('.study-col').toggle();
                $('.dt-buttons').children().first().text(function(i, text){
                      return text === "Hide study info" ? "Show study info" : "Hide study info";
                  });
                }},
             {extend: 'collection',
              text: 'Select company',
              autoClose: true,
              buttons: companies},
             {extend: 'collection',
              text: 'Select indication',
              autoClose: true,
              buttons: indications},
        ],
    });
});
function populateTable(data) {
    $('#development_table tbody').empty();
        $.each(data.dev_data, function (key, val) {
           let end_date = (val.date) ? val.date : "";
           let title = (val.title) ? val.title : "";
           $('#development_table').append(
                '<tr><td>' + val.phase+'</td>' +
                '<td><a href=\"' + val.url + '\"><img src=' + val.structure_url + ' title=\"' + val.name + '\" height=\"70\" width=\"70\"></a>' + '</td>' +
                '<td><a href=\"' + val.url + '\">' + val.name + '</a></td>' +
                '<td><a href=\"' + val.activity_url + '\">' + val.activity + '</a></td>' +
                '<td>' + val.company + '</td>' +
                '<td class=\"study-col\" style=\"display: none\">' + end_date + '</td>' +
                '<td class=\"study-col\" style=\"display: none\">' + title + '</td></tr>');
    });
}