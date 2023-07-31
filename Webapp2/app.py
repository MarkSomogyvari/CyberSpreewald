from shiny import *#App, render, ui
from main1 import model_channels, plot_stuff

app_ui = ui.page_fluid(
    ui.head_content(
        ui.row(
            ui.input_switch("w1", "Open weir", value = False),
        ),
        ui.row(
            ui.output_plot("pie"), #width="500px", height="1000px"
        ),
    ),
)

def server(input, output, session):
    @output
    @render.plot(alt="A histogram")
    
    def pie():


        sf, P_new, V = model_channels(W=input.w1())
        fig = plot_stuff(sf, P_new, V)

        return fig




app = App(app_ui, server, debug=True)