from dash import Dash, html, dash_table, dcc, Output, Input
import dash_bootstrap_components as dbc
import numpy as np
from math import sqrt
import pathlib
import scanpy as sc
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import analysis as an


# Initialise Dash app
app = Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    prevent_initial_callbacks="initial_duplicate"
)
app.title = "SEACells - Single-Cell Expression Explorer for SEACells"

# Suppress callback exceptions
app.config.suppress_callback_exceptions = True

# Path
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()

# Load data
# TODO: Load data from user upload

# Load anndata containing single-cell data
ad = sc.read(DATA_PATH.joinpath("anndata.h5ad"))
# Load A, the assignment matrix of cells to SEACells
A = np.load(DATA_PATH.joinpath("A.npy"))

# Process anndata and assignment matrix to identify triangles which define phenotypes
data = an.triangles(ad, A)


    
def description_card():
    """

    :return: A Div containing dashboard title & descriptions.
    """
    return html.Div(
        id="description-card",
        style={"width": "100%"},
        children=[
            html.H5("SEECells"),
            html.H3("Welcome to the SEACells Visualization Dashboard!"),
            html.Div(
                id="intro",
                children="Explore single cell data in the space of proximity to SEACells.",
            ),
        ],
    )

old_clickData = None
def generate_SEACells_description(clickData = None, triangle_no = 0):
    global old_clickData
    if clickData is None and old_clickData is None:
        tri_adj_matrix, tri_labels, tri_coords, strength = an.triangle_info(ad, A)
        triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength)
        triangles = []
        return triangle_fig
    if(clickData is None):
        clickData = old_clickData
    else:
        old_clickData = clickData

    clicked_seacell_no =  clickData['points'][0]['pointIndex']
    confirmed_triangles = data.confirmed_triangles
    triangles = [triangle for triangle in confirmed_triangles if triangle[0] == clicked_seacell_no or triangle[1] == clicked_seacell_no or triangle[2] == clicked_seacell_no]
    
    if len(triangles) == 0:
        tri_adj_matrix, tri_labels, tri_coords, strength = an.triangle_info(ad, A)
        triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength)
        return triangle_fig

    seacell_1 = triangles[triangle_no][0]
    seacell_2 = triangles[triangle_no][1]
    seacell_3 = triangles[triangle_no][2]

    tri_adj_matrix, tri_labels, tri_coords, strength = an.triangle_info(ad, A, seacell_1, seacell_2, seacell_3)

    triangle_fig  = triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength)
    return triangle_fig

def SEACells_description_card():
    """
    Construct a Div containing SEACells descriptions for a given SEACell. This is used to populate the
    SEACells description card. If no SEACell is selected, the card is empty.
    :param SEACell: (str) The name of the SEACell to describe. If None, the card is empty.
    :return: A Div containing SEACells description.
    """
    triangle_fig = generate_SEACells_description()

    return html.Div(
        id="SEACells-description-card",
        style={"height": "100%", "width": "100%"},
        children=[
            html.H1("SEACells Triangle Layout"),
            dcc.Dropdown(
                id='triangles_dropdown',
                value=1,
                placeholder="Select a triangle to view",
            ),
            dcc.Graph(id='SEACells-description-graph',
                style={"height": "100%", "width": "100%"},
                figure=triangle_fig
            ),
        ]
    )

@app.callback(
    Output('SEACells-description-graph', 'figure', allow_duplicate=True),
    Input('triangles_dropdown', 'value')
)
def update_triangles(value):
    print(f'Update to {value}')
    return generate_SEACells_description(triangle_no=value-1)

@app.callback(
    Output('triangles_dropdown', 'options'),
    [Input('SEACells-network-graph', 'clickData')])
def set_dropdown_options(clickData):
    if(clickData is None):
        return []
    clicked_seacell_no = clickData['points'][0]['pointIndex']
    list_of_triangles_for_each_seacell = data.list_of_triangles_for_each_seacell
    if clicked_seacell_no >= len(list_of_triangles_for_each_seacell):
        return []
    triangles = list_of_triangles_for_each_seacell[clicked_seacell_no]
    column_options = []
    for i, triangle in enumerate(triangles):
        column_options.append({'label': f'SEACell-{triangle[0]}, SEACell-{triangle[1]}, SEACell-{triangle[2]}', 'value': i+1})
    return column_options
'''
@app.callback(
    Output('triangles-dropdown', 'value'),
    [Input('triangles-dropdown', 'options')])
def set_cities_value(available_triangles):
    return available_triangles[0]['value']
'''



def triangle_layout(tri_adj_matrix, tri_labels, tri_coords, strength):

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
                title='Strength',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))

    node_strength = []
    node_text = []

    for i, node in enumerate(G.nodes()):
        if G.nodes[node]['label'].startswith("SEACell-"):
            continue
        node_strength.append(strength[i])
        node_text.append(f'This is {G.nodes[node]["label"]} <br># of strength: {strength[i]} ')

    node_trace.marker.color = node_strength
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
    :param selected_SEACell: (str) The name of the SEACell to highlight in the UMAP plot. If None, no SEACell is
                            highlighted. Otherwise, the SEACell is highlighted in red.
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
    :param selected_SEACell: (str) The name of the SEACell to highlight in the UMAP plot. If None, no SEACell is 
                                highlighted. Otherwise, the SEACell is highlighted in red.

    :return: A Div containing a UMAP plot of the single-cell data.
    """
    # Create dropdown options
    column_options = [{'label': i, 'value': i} for i in ad.obs.columns]
    if 'celltype' in ad.obs.columns:
        default = 'celltype'
    else:
        default = column_options[0]['value']

    # Generate UMAP plot
    cell_umap_fig = plot_sc_umap(default, selected_SEACell)
    return html.Div(
        id="sc-umap-card",
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
                style={"height": "100%", "width": "100%"}
            ),
        ],
        style={"height": "100%", "width": "100%"}
    )

@app.callback(
    Output('sc-umap-plot', 'figure'),
    Input('sc_umap_dropdown', 'value')
)
def update_sc_umap_plot(value):
    print(f'Colouring sc UMAP plot by {value}')
    return plot_sc_umap(value)


def plot_network_layout(adjacency_matrix, seacell_labels, UMAP_coords, clicked_point=None):
    """
    Plot a network layout of the SEACells anndata
    :param adjacency_matrix: The adjacency matrix of the SEACells anndata
    :param seacell_labels: The labels of the SEACells anndata
    :param UMAP_coords: The UMAP coordinates of the SEACells anndata
    :param clicked_point: The point clicked on the UMAP plot. If None, no point is highlighted.
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
                    marker=dict(
                        size=12.5,
                        color='green',
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
# (3) the SEACells description (to show the description of the selected SEACell on click)
@app.callback(
    [
        Output('SEACells-network-graph', 'figure'),
        Output('sc-umap-plot', 'figure', allow_duplicate=True),
        Output('SEACells-description-graph', 'figure')
        #Output('triangles_dropdown', 'options', 'value')
    ],
    Input('SEACells-network-graph', 'clickData')
)
def update_network_layout_and_highlight_cells_on_sc_umap(clickData):
    print("ENTERED_CALLBACK_MESSAGE")
    """
    Update the network layout to highlight the selected SEACell on click and the UMAP plot to highlight the cells
    assigned to the selected SEACell on click.
    :param clickData: The data of the point clicked on the UMAP plot.
    """
    # Extract the selected SEACell from the clickData
    if clickData is not None:
        selected_SEACell = 'SEACell-' + str(clickData['points'][0]['pointIndex'])
    else:
        selected_SEACell = None

    print(f'Updating network layout and highlighting cells assigned to SEACell {selected_SEACell} on click.')
    return plot_network_layout(data.adjacency_matrix, data.seacell_labels, data.UMAP_coords, clickData), \
           plot_sc_umap(None, selected_SEACell), \
           generate_SEACells_description(clickData)
           #SEACells_description_card(clickData)


# calls plot_network_layout to construct the figure
def SEACells_network_card():
    """
    Generate a div containing a network layout of the SEACells.

    """
    SEACells_network_fig = plot_network_layout(data.adjacency_matrix, data.seacell_labels, data.UMAP_coords)
    return html.Div(
        id="SEACells-network-card",
        style={"height": "100%", "width": "100%"},
        children=[
            html.H1("SEACells Network Layout"),
            dcc.Graph(id='SEACells-network-graph',
                      style={"height": "100%", "width": "100%"},
                      figure=SEACells_network_fig)
        ],
    )

app.layout = html.Div(
    id="app-container",
    children=[
        # Left column
        html.Div(
            id="left-column",
            className="five columns",
            children=[
                description_card(),
                sc_umap_card(),
            ]
        ),
        # Right column
        html.Div(
            id="right-column",
            className="seven columns",
            children=[
                SEACells_network_card(),
                SEACells_description_card(),
            ],
        ),
    ],
)


# Run the server
if __name__ == "__main__":
    app.run_server(debug=True)
