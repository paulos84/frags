$(document).ready(function(){

    $("#id_property_choice").change(function () {
      var choice = $("#id_property_choice").val();
      $(".chem_data_plot").hide();
      $("#" + choice).show();
    });

    $("#show-odor").click(function () {
      $(this).text(function(i, text){
          return text === "Synonyms" ? "Odor types" : "Synonyms";
      })
      $('.odor-terms').toggle();
      $('.compound-synonyms').toggle();
    });
});
