import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
from dash import Dash, html, dash_table, dcc, Output, Input, State, clientside_callback
import dash_bootstrap_components as dbc
import numpy as np
from math import sqrt
import pathlib
import scanpy as sc
import pandas as pd
import analysis as an
import dash
import SEACells as SEACells

import dash_uploader as du
from time import sleep
from flask import Response
import webbrowser
from numpy import linalg

import time

from worker import seacells_algorithm
from email_sender import send_email

from redis import Redis
from rq import Queue


# Initialise Dash app
app = Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    prevent_initial_callbacks="initial_duplicate",
    use_pages=True,
)
app.title = "SEACells - Single-Cell Expression Explorer for SEACells"

# Suppress callback exceptions
app.config.suppress_callback_exceptions = True

# Path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()
du.configure_upload(app, r"upload")




def description_card():
    """
    :return: A Div containing dashboard title & descriptions.
    """
    style = { 
                "height": "auto",
                "padding": "auto",
                "margin": "auto",
                "width": "30%",
                'text-align': 'center',
            }
    return html.Div(
        id="description-card",
        style={"width": "100%"},
        children=[
            html.Div(
    style={'display': 'flex', 'justifyContent': 'center', 'alignItems': 'center'},  # This centers the content and aligns items horizontally
    children=[
        html.A(
            html.Img(
                src=app.get_asset_url('home.jpg'),
                style={'height': '50px', 'width': '50px', 'marginRight': '10px'}  # marginRight adds space between the image and header
            ),
            href='/',
            style={'alignSelf': 'center'}  # This ensures the image is centered vertically within the flex container
        ),
        dcc.Link(
            html.H5(
                "SEACells",
                style={'textAlign': 'center', 'font-size': '4.2rem', 'cursor': 'pointer', 'marginBottom': '0px'}
            ),
            href="/",
            style={'text-decoration': 'none', 'display': 'block'}  # display block is added to ensure proper alignment within the flex container
        )
    ]
),
            html.H3("Welcome to the SEACells Visualization Dashboard!", style={'textAlign': 'center', 'marginTop': '0px'}),
        ],
    )

@app.callback(
    [Output("email-input", "valid"), Output("email-input", "invalid")],
    [Input("email-input", "value")],
)
def check_validity(text):
    if text:
        is_gmail = text.endswith("@gmail.com")
        return is_gmail, not is_gmail
    return False, False
            
    

@du.callback(
    output=[Output("display-val-container", "hidden"), Output("filename", "value")],
    id="cell-anndata",
)
def single_cell_callback(status: du.UploadStatus):
    fname = status.latest_file
    return False, str(fname)


@app.callback(
    Output(component_id='main-holder', component_property='hidden'),
    Output(component_id='msg', component_property='children'),
    Output(component_id='msg', component_property='hidden'),
    Output(component_id='display-val', component_property='disabled'), # Can change later! temporary solution
    Input('display-val','n_clicks'),
    State('filename', 'value'),
    State('num-seacells', 'value'),
    State('email-input', 'value'),
    prevent_initial_call=True)
def update(n_clicks, fname, num_seacells, email_address):
    if (n_clicks >= 1):
        if (not email_address.endswith('@gmail.com')):
            return False, [html.P(f"Correct the email address your provided: {email_address}", style={'textAlign': 'center'})], False, None
        
        ad = sc.read(fname)
        no_cells = sc.AnnData(ad.X).shape[0]
        minutes = 2
        return True, [html.P(f"Your dataset {fname.split('/')[-1]} is successfully submitted. Your dataset contains {no_cells} many cells, and you want {num_seacells} many SEACells.", style={'textAlign': 'center'}), html.P(f"The SEACELLS algorithm will take approximately {minutes} minutes. The results will be emailed to your address.", style={'textAlign': 'center'})], False, 'disabled'
    else:
        return False, [], True, None

@app.callback(
    Output('location', 'href'),
    Input('display-val','n_clicks'),
    State('filename', 'value'),
    State('num-seacells', 'value'),
    State('email-input', 'value'),
    prevent_initial_call=True)
def update_visibility(n_clicks, fname, num_seacells, email_address):
    if (n_clicks >= 1):
        if (not email_address.endswith('@gmail.com')):
            return None

        r = Redis(host='localhost', port=6379)
        q = Queue(connection=r)

        job = q.enqueue(seacells_algorithm, result_ttl = 259200,
                        job_timeout=1500, args=[fname, num_seacells])

        while(job.return_value() is None):
            print('waiting')
            time.sleep(2)

        URL = 'http://127.0.0.1:8050/' + job.id
        subject = "Your SEACells processing is complete."
        message = "Please visit " + URL +  " to view the results."
        send_email(subject, message, email_address)
        return URL
    else:
        return None

app.layout = html.Div(
    id="app-container",
    children=[
        dcc.Location(id='location'),
        html.Div(children=[description_card()]),
        dash.page_container,
    ],
)


# Run the server
if __name__ == "__main__":
    app.run_server(debug=False)
    webbrowser.open('http://127.0.0.1:8050/')