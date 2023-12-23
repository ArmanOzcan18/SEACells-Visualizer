import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")
from dash import Dash, html, dash_table, dcc, Output, Input, State, callback
import dash_bootstrap_components as dbc
import numpy as np
from math import sqrt
import pathlib
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import analysis as an
import dash
import SEACells as SEACells

import dash_uploader as du
from time import sleep
from flask import Response
from numpy import linalg

from celery import Celery
import pickle

from worker import seacells_algorithm

from email_sender import send_email

from redis import Redis
from rq import Queue

dash.register_page(__name__, path_template="/<report_id>")

# Suppress callback exceptions
# app.config.suppress_callback_exceptions = True

# Path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()

# du.configure_upload(app, r"upload")

# ad = pickle.load(open('pickle/ad.pickle', 'rb'))
# SEACell_ad = pickle.load(open('pickle/SEACELL_ad.pickle', 'rb'))
# A = pickle.load(open('pickle/A.pickle', 'rb'))


r = Redis(host='localhost', port=6379)
q = Queue(connection=r)

ad, SEACell_ad, A, data = None, None, None, None

#data = pickle.load(open('pickle/data.pickle', 'rb'))



def generate_SEACells_description(clickData = None, triangle_no = 0, feature_no = 0):
    """
    Construct a Div containing SEACells triangle layout for a given SEACell.
    :param clickData: The data of the point clicked on the UMAP plot.
    :param triangle_no: The triangle number to display. If not specified, the first triangle is displayed.
    :return: A Div containing SEACells description.
    """
    tri_adj_matrix, tri_labels, tri_coords, strength = an.triangle_info(ad, A)

    if feature_no == 0:
        triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength = strength)
    else:
        triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, gene = ad.var_names[feature_no-1])
    
    if clickData is None:
        return triangle_fig
    
    clicked_seacell_no =  clickData['points'][0]['pointIndex']
    confirmed_triangles = data.confirmed_triangles
    triangles = [triangle for triangle in confirmed_triangles if triangle[0] == clicked_seacell_no or triangle[1] == clicked_seacell_no or triangle[2] == clicked_seacell_no]
    
    if (len(triangles) == 0):
        return triangle_fig

    seacell_1 = triangles[triangle_no][0]
    seacell_2 = triangles[triangle_no][1]
    seacell_3 = triangles[triangle_no][2]

    tri_adj_matrix, tri_labels, tri_coords, strength = an.triangle_info(ad, A, seacell_1, seacell_2, seacell_3)

    if feature_no == 0:
        triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength = strength)
    else:
        triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, gene = ad.var_names[feature_no-1])
    
    return triangle_fig

def SEACells_description_card():
    """
    Construct a Div containing SEACells descriptions for a given SEACell.
    :return: A Div containing SEACells description.
    """
    triangle_fig = generate_SEACells_description()

    feature_options = [{'label': "Assignment Strength", 'value': 0}]
    for i, feature in enumerate(ad.var_names):
        feature_options.append({'label': feature, 'value': i+1})

    return html.Div(
        id="SEACells-description-card",
        style={"height": "100%", "width": "100%"},
        children=[
            html.H5("SEACells Triangle Layout"),
            dcc.Dropdown(
                id='triangles_dropdown',
                value=0,
                placeholder="Select a triangle to view",
            ),
            dcc.Dropdown(
                id='features_dropdown',
                value=0,
                options=feature_options,    
                placeholder="Select a feature to view",
            ),
            dcc.Graph(id='SEACells-description-graph',
                style={"height": "100%", "width": "100%"},
                figure=triangle_fig
            ),
        ]
    )

@callback(
    Output('SEACells-description-graph', 'figure', allow_duplicate=True),
    [Input('triangles_dropdown', 'value'),
     Input('features_dropdown', 'value'),
     Input('SEACells-network-graph', 'clickData')], prevent_initial_call=True
)
def update_triangles(triangle_no, feature_no, clickData):
    """
    Update the SEACells description to show the selected triangle and feature.
    :param triangle_no: The triangle to display.
    :param feature_no: The feature of the triangle to display
    :param clickData: The data of the point clicked on the UMAP plot.
    :return: The updated SEACells description.
    """
    if clickData is None:
        return generate_SEACells_description(feature_no=feature_no)
    
    if triangle_no is None:
        print('ENTERED TRIANGLE NO IS NONE')
        return generate_SEACells_description(clickData = clickData, triangle_no=0, feature_no=feature_no)
    
    clicked_seacell_no = clickData['points'][0]['pointIndex']

    if len(data.list_of_triangles_for_each_seacell) <= clicked_seacell_no:
        print('How can this happen?')
        return generate_SEACells_description(clickData = clickData, triangle_no=0, feature_no=feature_no)
    
    triangles = data.list_of_triangles_for_each_seacell[clicked_seacell_no]

    if len(triangles) <= triangle_no:
        return generate_SEACells_description(clickData = clickData, triangle_no=0, feature_no=feature_no)
    
    print('Updated the triangle')

    return generate_SEACells_description(clickData = clickData, triangle_no=triangle_no, feature_no=feature_no)





@callback(
    Output('triangles_dropdown', 'options'),
    [Input('SEACells-network-graph', 'clickData')])
def set_dropdown_options(clickData):
    """
    Set the options of the triangles dropdown to the triangles of the selected SEACell.
    :param clickData: The data of the point clicked on the UMAP plot.
    :return: The options of the triangles dropdown.
    """
    if(clickData is None):
        return []
    clicked_seacell_no = clickData['points'][0]['pointIndex']

    if len(data.list_of_triangles_for_each_seacell) <= clicked_seacell_no:
        print("ENTERED!!!") # I don't know why this is happening sometimes!
        return []

    triangles = data.list_of_triangles_for_each_seacell[clicked_seacell_no]
    column_options = []
    for i, triangle in enumerate(triangles):
        column_options.append({'label': f'SEACell-{triangle[0]}, SEACell-{triangle[1]}, SEACell-{triangle[2]}', 'value': i})

    return column_options
'''
@app.callback(
    Output('triangles-dropdown', 'value'),
    [Input('triangles-dropdown', 'options')])
def set_cities_value(available_triangles):
    return available_triangles[0]['value']
'''

def triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength = None, gene = None):
    """
    Plot a triangle of SEACells
    :param adjacency_matrix: The adjacency matrix of SEACells on the vertices of the triangle
    :param seacell_labels: The labels of the cells
    :param UMAP_coords: The UMAP coordinates of the cells
    :param strength: The strength of each cell in that triangle; in other words, how assigned that cell is to the triangle
    :return: A triangle figure of the SEACells
    """
    G = an.adjacency_matrix_to_graph(tri_adj_matrix)
    G = an.annotate_nodes(tri_labels, G)
    G = an.add_coordinates(tri_coords, G)

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['x'], G.nodes[edge[0]]['y']
        x1, y1 = G.nodes[edge[1]]['x'], G.nodes[edge[1]]['y']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    node_label = []

    if strength is not None:
        feature_title = 'Assignment Strength'
        node_color = 'Hot'
    else:
        feature_title = gene.capitalize() + ' Expression'
        node_color = 'Earth'

    for node in G.nodes():
        x, y = G.nodes[node]['x'], G.nodes[node]['y']
        node_x.append(x)
        node_y.append(y)
        node_label.append(G.nodes[node]['label'])

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale=node_color,
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title=feature_title,
                xanchor='left',
                titleside='right'
            ),
            line_width=2))


    node_text = []
    node_strength = []
    gene_expression = []

    if strength is None:
        cell_gene_series = ad.to_df().loc[:, gene]
        SEACell_gene_series = SEACell_ad.to_df().loc[:, gene]

    for i, node in enumerate(G.nodes()):
        name = G.nodes[node]["label"]

        if strength is not None:
            node_strength.append(strength[i])
            
            if(name.startswith("SEACell")):
                node_text.append(name)
            else:
                node_text.append(f'{name} of type {ad.obs["celltype"][name]} of strength {strength[i]}')
        else:
            if(name == 'SEACell'):
                gene_expression.append(0)
                node_text.append(name)
            elif(name.startswith("SEACell")):
                gene_expression.append( SEACell_gene_series[name] )
                node_text.append(f'{name} of {gene} expression { SEACell_gene_series[name] }')
            else:
                gene_expression.append( cell_gene_series[name] )
                node_text.append(f'{name} of type {ad.obs["celltype"][name]} of  {gene} expression: {cell_gene_series[name]}')

    if strength is not None:
        node_trace.marker.color = node_strength
    else:    
        node_trace.marker.color = gene_expression
    
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        annotations=[dict(
                            text="",
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002)],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    
    fig.add_trace(go.Scatter(
        hoverinfo='skip',
        x=[0, 1, -1],
        y=[sqrt(3)+0.05, -0.125, -0.125],
        mode="text",
        name="Lines, Markers and Text",
        text=[tri_labels[0], tri_labels[1], tri_labels[2]],
        textposition="top center")
        )
    return fig




def plot_sc_umap(column, selected_SEACell=None):
    """
    Plot a UMAP plot of the single-cell data.
    :param column: (str) The column to color the UMAP plot by.
    :param selected_SEACell: (str) The name of the SEACell to highlight in the UMAP plot. If None, no SEACell is highlighted.
    :return: A UMAP plot of the single-cell data.
    """
    if column is None:
        if 'celltype' in ad.obs.columns:
            column = 'celltype'
        else:
            column = ad.obs.columns[0]

    print(f'Coloring by {column} when updating UMAP plot.')
    cell_umap_fig = px.scatter(ad.obsm['X_umap'],
                               x=0,
                               y=1,
                               color=ad.obs[column],
                               hover_name=ad.obs.index,
                               labels={'color': column},
                             )
    # Update all traces to have the same marker size and opacity
    cell_umap_fig.update_traces(marker={'size': 2})

    if selected_SEACell is not None:
        print(f'Highlighting {selected_SEACell} in UMAP plot.')
        # Get numeric indices of cells in the selected SEACell
        selected_cells = np.where(ad.obs['SEACell'] == selected_SEACell)[0]
        # Add a trace outlining the selected cells in black
        cell_umap_fig.add_trace(go.Scatter(x=ad.obsm['X_umap'][selected_cells, 0],
                                             y=ad.obsm['X_umap'][selected_cells, 1],
                                             mode='markers',
                                             marker=dict(size=5,
                                                         line=dict(width=2,
                                                                    color='darkslategray')
                                                         ),
                                             showlegend=False,
                                             hoverinfo='none'
                                             ))

    return cell_umap_fig

def sc_umap_card(selected_SEACell=None):
    """
    Generate a div containing a UMAP plot of the single-cell data.
    :param selected_SEACell: (str) The name of the SEACell to highlight in the UMAP plot. If None, no SEACell is highlighted.
    :return: A Div containing a UMAP plot of the single-cell data.
    """

    # Create dropdown options
    column_options = [{'label': i, 'value': i} for i in ad.obs.columns]
    if 'celltype' in ad.obs.columns:
        default = 'celltype'
    else:
        default = column_options[0]['value']

    # Generate UMAP plot
    print('Generating UMAP plot!!!')
    cell_umap_fig = plot_sc_umap(default, selected_SEACell)

    return html.Div(
        id="sc-umap-card",
        style={"width": "100%"},
        children=[
            html.H5("Single-cell UMAP"),
            dcc.Dropdown(
                id='sc_umap_dropdown',
                options=column_options,
                value=default
            ),
            dcc.Graph(
                id="sc-umap-plot",
                figure=cell_umap_fig,
            )
        ]
    )

@callback(
    Output('sc-umap-plot', 'figure'),
    Input('sc_umap_dropdown', 'value'), prevent_initial_call=True)
def update_sc_umap_plot(value):
    """ 
    Update the UMAP plot to color by the selected column.
    :param value: The column to color the UMAP plot by.
    :return: The updated UMAP plot.
    """
    print(f'Colouring sc UMAP plot by {value}')
    return plot_sc_umap(value)



def plot_network_layout(adjacency_matrix, seacell_labels, UMAP_coords, clicked_point=None):
    """
    Plot a network layout of the SEACells
    :param adjacency_matrix: The adjacency matrix of the SEACells network graph
    :param seacell_labels: The labels of the SEACells
    :param UMAP_coords: The UMAP coordinates of the SEACells
    :param clicked_point: The point clicked on the UMAP plot. If None, no point is highlighted.
    :return: A network layout of the SEACells
    """
    G = an.adjacency_matrix_to_graph(adjacency_matrix)
    G = an.annotate_nodes(seacell_labels, G)
    G = an.add_coordinates(UMAP_coords, G)

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['x'], G.nodes[edge[0]]['y']
        x1, y1 = G.nodes[edge[1]]['x'], G.nodes[edge[1]]['y']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    node_label = []

    for node in G.nodes():
        x, y = G.nodes[node]['x'], G.nodes[node]['y']
        node_x.append(x)
        node_y.append(y)
        node_label.append(G.nodes[node]['label'])

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append(f'SEACell: {node_label[node]}<br># of connections: {len(adjacencies[1])}')

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        annotations=[dict(
                            text="",
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002)],
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    if clicked_point is not None:
        print('Updating network layout with highlighted SEACell')
        clicked_seacell_no = clicked_point['points'][0]['pointIndex']
        confirmed_triangles = data.confirmed_triangles
       
        triangles = [triangle for triangle in confirmed_triangles if triangle[0] == clicked_seacell_no or triangle[1] == clicked_seacell_no or triangle[2] == clicked_seacell_no]
        
        neighbor_labels = set(np.array(triangles).flatten())
        if(clicked_seacell_no in neighbor_labels):
            neighbor_labels.remove(clicked_seacell_no)       

        fig.add_trace(
            go.Scatter(
                x=[clicked_point['points'][0]['x']],
                y=[clicked_point['points'][0]['y']],
                mode='markers',
                hoverinfo='none',
                marker=dict(
                    size=20,
                    color='red',
                    line=dict(
                        width=2,
                        color='black'
                    )
                )
            )
        )
        
        for neighbor_label in neighbor_labels:
            fig.add_trace(
                go.Scatter(
                    x=[G.nodes[neighbor_label]['x']],
                    y=[G.nodes[neighbor_label]['y']],
                    mode='markers',
                    hoverinfo='none',
                    marker=dict(
                        size=12.5,
                        color='orange',
                        line=dict(
                            width=2,
                            color='black'
                        )
                    )
                )
            )

    return fig

# Use callbacks to update:
# (1) the network layout (to highlight the selected SEACell on click)
# (2) the UMAP plot (to highlight the cells assigned to the selected SEACell on click).
@callback(
    [
        Output('SEACells-network-graph', 'figure'),
        Output('sc-umap-plot', 'figure', allow_duplicate=True)
    ],
    Input('SEACells-network-graph', 'clickData'), prevent_initial_call=True
)
def update_network_layout_and_highlight_cells_on_sc_umap(clickData):
    """
    Update the network layout to highlight the selected SEACell on click and the UMAP plot to highlight the cells
    assigned to the selected SEACell on click.
    :param clickData: The data of the point clicked on the UMAP plot.
    :return: The updated network layout and UMAP plot.
    """
    # Extract the selected SEACell from the clickData
    if clickData is not None:
        selected_SEACell = 'SEACell-' + str(clickData['points'][0]['pointIndex'])
    else:
        selected_SEACell = None

    print(f'Updating network layout and highlighting cells assigned to SEACell {selected_SEACell} on click.')
    return plot_network_layout(data.adjacency_matrix, data.seacell_labels, data.UMAP_coords, clickData), \
           plot_sc_umap(None, selected_SEACell), \
           

# calls plot_network_layout to construct the figure
def SEACells_network_card():
    """
    Generate a div containing a network layout of the SEACells.
    :return: A Div containing a network layout of the SEACells.
    """

    SEACells_network_fig = plot_network_layout(data.adjacency_matrix, data.seacell_labels, data.UMAP_coords)

    return html.Div(
        id="SEACells-network-card",
        style={"height": "100%", "width": "100%"},
        children=[
            html.H5("SEACells Network Layout"),
            dcc.Graph(id='SEACells-network-graph',
                      style={"height": "100%", "width": "100%"},
                      figure=SEACells_network_fig)
        ],
    )



def boxplot_card():
    """
    Generate a div containing a boxplot of some statistics.
    :return: A Div containing the boxplot.
    """
    list = ["purity", "distplot"]
    # list = ["purity", "compactness", "separation", "distplot"]
    # Create dropdown options
    column_options = [{'label': i, 'value': i} for i in list]
    default = "purity"
    # Generate UMAP plot
    fig = plot_boxplots(default)

    return html.Div(
        id="boxplot-card",
        children=[
            html.H5("SEACell Boxplots"),
            dcc.Dropdown(
                id='boxplot_dropdown',
                options=column_options,
                value=default
            ),
            dcc.Graph(
                id="boxplot",
                figure = fig,
                style={"height": "100%", "width": "100%"}
            ),
        ],
        style={"height": "100%", "width": "100%"}
    )

@callback(
    Output('boxplot', 'figure'),
    Input('boxplot_dropdown', 'value')
)
def update_boxplot(value):
    """
    Update the boxplot to show the selected statistic.
    :param value: The statistic to show.
    :return: The updated boxplot.
    """
    print(f'Changing boxplot to {value}')
    return plot_boxplots(value)
    
def plot_boxplots(version):
    """
    Plot a boxplot of some statistics.
    :param version: The version of the boxplot to plot.
    :return: A boxplot of some statistics.
    """
    if(version == "purity"):
        df = SEACells.evaluate.compute_celltype_purity(ad, 'celltype')
        fig = px.box(df, y="celltype_purity", title = "Celltype Purity")
        return fig
    elif(version == "compactness"):
        df = SEACells.evaluate.compactness(ad, 'X_pca')
        fig = px.box(df, y="compactness", title = "Compactness of SEACells")
        return fig
    elif(version == "separation"):
        df = SEACells.evaluate.separation(ad, 'X_pca',nth_nbr=1)
        fig = px.box(df, y="separation", title = "Separation of SEACells")
        return fig
    elif(version == "distplot"):
        values = (A.T > 0.1).sum(axis=1)
        df = pd.DataFrame({'number of SEACells': values})
        fig = px.histogram(df, x="number of SEACells", title = "Non-trivial (>0.1) Assignments per Cell", nbins=7)
        return fig
    else:
        print("NOOOO!")
        return None






def layout(report_id=None):

    global ad, SEACell_ad, A, data
    job = q.fetch_job(report_id)
    ad, SEACell_ad, A, data = job.return_value()

    return html.Div(
    id="app-container",
    style={},
    children=[
        html.Div(id="left-column",
            className="six columns",
            children=[sc_umap_card(),boxplot_card()],
        ),
        # Right column
        html.Div(id="right-column",
            style={"background-color": "#f9f9f9",},
            className="six columns",
            children=[SEACells_network_card(), SEACells_description_card()],
        ),
    ],
)