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
   $('.your-class').slick({
      infinite: true,
      slidesToShow: 6,
      slidesToScroll: 6
    });
});
function initializeTable(data) {
    var table_rows = []
    $.each(data, function (key, val) {
       let end_date = (val.date) ? val.date : "";
       let title = (val.title) ? val.title : "";
        table_rows.push([val.phase,
        '<a href=\"' + val.url + '\"><img src=' + val.structure_url + ' title=\"' + val.name + '\" height=\"70\" width=\"70\"></a>',
                '<a href=\"' + val.url + '\">' + val.name + '</a>',
                '<a href=\"' + val.activity_url + '\">' + val.activity + '</a>',
                val.company,
                end_date,
                title])
    });
    $('#development_table').dataTable({
        paging: true,
        bAutoWidth: false ,
        dom: 'Bfrtip',
        buttons: [
             {extend: 'collection',
              text: 'Select company',
              autoClose: true,
              buttons: companies},
             {extend: 'collection',
              text: 'Select category',
              autoClose: true,
              buttons: categories},
             {extend: 'collection',
              text: 'Select indication',
              autoClose: true,
              buttons: indications},
              ],
          data: table_rows,
          columns: [
            { title: "Phase" },
            { title: "Structure" },
            { title: "Compound" },
            { title: "Indication/Mechanism" },
            { title: "Company" },
            { title: "Ends" },
            { title: "Study Title" }
        ],
        columnDefs: [
            { className: "hide_column", targets: [ 5 ] },
            { className: "hide_column", targets: [ 6 ] },
            {bSortable: false, targets: [ 1 ]  },
            {bSortable: false, targets: [ 6 ]  }
        ],
        order: [[ 0, "desc" ]],
        aLengthMen: [[10, 25, 50, -1], [10, 25, 50, "All"]]
    });
}
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
            '<td class=\"study-col\" >' + end_date + '</td>' +
            '<td class=\"study-col\" >' + title + '</td></tr>');
});
}