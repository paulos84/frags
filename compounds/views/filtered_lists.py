"""    class BookListView(generic.ListView):
        model = Book
        context_object_name = 'my_book_list'   # your own name for the list as a template variable
        queryset = Book.objects.filter(title__icontains='war')[:5] # Get 5 books containing the title war
        template_name = 'books/my_arbitrary_template_name_list.html'  # Specify your own template name/location

     Override the get_queryset() method to change the list of records returned. This is more flexible than just setting
     the queryset attribute as in the preceding code fragment

    class BookListView(generic.ListView):
    model = Book
        def get_queryset(self):
            return Book.objects.filter(title__icontains='war')[:5] # Get 5 books containing the title war

     Override get_context_data() in order to pass additional context variables to the template. The code below shows how
     to add a variable named "some_data" to the context (it would then be available as a template variable).

     """