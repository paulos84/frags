$(document).ready(function() {

$( "#clickme" ).click(function() {
  $( ".breakword" ).hide( "slow", function() {
    alert( "Animation complete." );
  });
});

});