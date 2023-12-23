import dash
import dash_daq as daq
import dash_uploader as du
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html

dash.register_page(__name__, path='/')

style = { 
                "height": "auto",
                "padding": "auto",
                "margin": "auto",
                "width": "30%",
                'text-align': 'center',
            }

layout = html.Div([
    html.Div(
        id="description-card",
        style={"width": "100%"},
        children=[
            html.Div(id='main-holder', hidden=False, children=[
                html.P("Explore single cell data in the space of proximity to SEACells.", style={'textAlign': 'center'}), 
                html.Div(id='upload', children=[
                    html.Div(className="twelve columns",
                        children=[
                            du.Upload(
                                id='cell-anndata',
                                text='Upload Cell Anndata',
                                # filetypes=['h5ad'],
                                default_style = style
                            ),
                            html.Div(id='display-val-container', hidden = True, style={'text-align': 'center'},
                                children=[
                                    html.Div(
                                        id='seacell_container',
                                        children=[
                                            html.Div(
                                                [
                                                    daq.NumericInput(
                                                        id='num-seacells',
                                                        style={'textAlign': 'center', 'marginBottom': '5px',  'marginRight': '0px'},
                                                        label='Number of SEACells',
                                                        min=1,
                                                        max=150,
                                                        labelPosition='top',
                                                        value=90,
                                                    ),
                                                    html.Div([dbc.Button(
                                                        "Tip",
                                                        id="tooltip-target-right",
                                                        n_clicks=0,
                                                    ),
                                                    dbc.Tooltip(
                                                        "The number of SEACells should be about 1/75 of the number of cells in the dataset.",
                                                        target="tooltip-target-right",
                                                        placement="right",
                                                    )], style={'marginTop': '25px', 'marginLeft': '0px'},
                                                    ),
                                                    html.Div(
                                                        [
                                                            dbc.Label("Email"),
                                                            dbc.Input(id="email-input", type="email",
                                                                    value="", style={'width': '300px'}),
                                                            dbc.FormText("We only accept gmail..."),
                                                            dbc.FormFeedback("That looks like a gmail address :-)", type="valid"),
                                                            dbc.FormFeedback(
                                                                "Sorry, we only accept gmail for some reason...",
                                                                type="invalid",
                                                            ),
                                                        ],
                                                        style={'marginLeft': '20px', 'marginTop' : '20px'},  # Space between the input fields
                                                    ),
                                                ],
                                                style={
                                                    'display': 'flex',
                                                    'justifyContent': 'center',
                                                    'alignItems': 'center'
                                                }
                                            )
                                        ],
                                        style={'textAlign': 'center'}
                                    ),
                                    dcc.Input(
                                        id='filename',
                                        type='hidden',
                                        value='',
                                    ),
                                    html.Button('Submit', id='display-val', n_clicks=0, style={"margin-top": "1%",}),
                                    ]
                            ),
                        ]
                    ),
                ]),
            ]),
            html.Div(id='msg', hidden=True, style={'text-align': 'center'}, children=[]),
        ],
    )
])