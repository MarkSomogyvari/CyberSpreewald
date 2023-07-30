from shiny import *#App, render, ui
from main1 import model_channels, plot_stuff

app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.input_switch("w1", "Open weir", value = False),
        ),
        ui.panel_main(
            ui.output_plot("pie"), width = '50%', height = '100%'
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