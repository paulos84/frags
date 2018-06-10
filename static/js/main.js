{
      var cas_number = $(this).val();
      $.ajax({
        url: '{% url "process_cas" %}',
        data: {
          'cas_number': cas_number
        },
        dataType: 'json',
        success: function (data) {
        alert(data)
          <!--if (data.is_taken) {-->
            <!--alert("A user with this username already exists.");-->
          }
        }
      });
    }