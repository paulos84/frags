from django.shortcuts import render
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import components

def example(request):
    plot = figure()
    plot.circle([1,2], [3,4])

    script, div = components(plot, CDN)

    return render(request, "simple_chart.html", {"the_script": script, "the_div": div})

