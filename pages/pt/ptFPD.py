import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import dash_table
from dash_table.Format import Format, Group, Scheme, Symbol

import dash_daq as daq
import dash_cytoscape as cyto
import sys
import os
current_dir = os.path.dirname(__file__)
sys.path.append(os.path.dirname(current_dir))
from pages.pt.ptUtils import Header, make_dash_table
import pandas as pd
import pathlib

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../../data").resolve()



def create_layout(app):
    return html.Div(
        [

            Header(app),
            html.Div(id='page-language', lang='pt-br', style={'display': 'none'}),
        
            html.Div([
                # The native browser confirm dialog
                dcc.ConfirmDialog(
                    id='confirm-error-fpd',
                        message='',
                ), 
                html.Div(id='confirm-fpd-success')
            ]),
        
            html.Div(
            [
            

                    html.H4(["Fluxo de Potência na Distribuição"], className="subtitle"),
                    html.Br(),
                    html.Div(
                    [
                        html.H6('Tolerância',style={'display':'inline-block','margin-right':20}),
                        html.Div(
                            [
                                dcc.Input(
                                    id="input_tol_FPD",
                                    type="number",
                                    value="1e-5",
                                    max = 1,
                                    min = 0.00001,
                                    step = 0.00001
                                )
                            ],
                        ),
                    ],className="six columns",
                    ),

                    html.Div(
                    [
                        html.H6('Limite Iterações',style={'display':'inline-block','margin-right':20}),
                        html.Div(
                            [
                                dcc.Input(
                                    id="input_iteracoes_FPD",
                                    type="number",
                                    value="50"
                                )
                            ],
                            
                        ),
                    ],className="six columns",
                    ),
                    
                    html.Div(
                    [
                        html.H6('Barra de Referência',style={'display':'inline-block','margin-right':20}),
                        html.Div(
                            [
                                dcc.Input(
                                    id="input_refbus_FPD",
                                    type="number",
                                    value="1",
                                )
                                # dcc.Dropdown(id='dynamic-dropdown', options=[], placeholder='Select an option'),
                                # html.Div(id='input_source_bus_FPD')
                            ],
                        ),
                    ],className="six columns",
                    ),
                    
                    html.Div(
                    [
                        html.H6('Barra final do ramo principal',style={'display':'inline-block','margin-right':20}),
                        html.Div(
                            [
                                dcc.Input(
                                    id="input_final_bus_FPD",
                                    type="number",
                                    value="1",
                                )
                            ],
                        ),
                    ],className="six columns",
                    ),


                    html.Div(
                    [
                        html.H6('Modulo Tensão Inicial (pu)',style={'display':'inline-block','margin-right':20}),
                        html.Div(
                            [
                                dcc.Input(
                                    id="input_v_FPD",
                                    type="number",
                                    # placeholder="1e-4",
                                    value="1",
                                    min = 0,
                                    step = 0.00001
                                )
                            ],
                        ),
                    ],className="six columns",
                    ),

                    html.Div(
                    [
                        html.H6('Fase Tensão Inicial (graus)',style={'display':'inline-block','margin-right':20}),
                        html.Div(
                            [
                                dcc.Input(
                                    id="input_ang_FPD",
                                    type="number",
                                    # placeholder="1e-4",
                                    value="0",
                                    max = 360,
                                    min = 0,
                                    step = 0.00001
                                )
                            ],
                        ),
                    ],className="six columns",
                    ),

                    html.Br(),
                    
                    html.Div(
                        [
                            html.H6(
                                ["Topologia"], className="subtitle padded"
                            ),
                            dcc.Upload(id='topology_FPD', disabled=True,
                                        children = html.Div([
                                                            html.A('Modifique os parametros no Editor de Topologia')
                                                            ]),),
                            dash_table.DataTable(id='topology_table_FPD', editable=False,row_deletable=False,page_action='none',
                            style_table={'height': '200px', 'overflowY': 'auto'},
                            style_cell={'textAlign': 'center'}),
                        ],
                        className="six columns",
                    ),
                        
                    html.Div(
                        [
                            html.H6(
                                ["Dados das Cargas"],
                                className="subtitle padded",
                            ),
                            dcc.Upload(id='load_FPD',
                                        children = html.Div([
                                                            html.A('Arraste ou Selecione o Arquivo')
                                                            ]),),
                            dash_table.DataTable(id='load_table_FPD',
                            editable=True,row_deletable= False,page_action='none',
                            style_table={'height': '200px', 'overflowY': 'auto'},
                            style_cell={'textAlign': 'center'}, selected_rows=list()),
                        ],
                        className="six columns",
                    ),
                    
                    
                    
                    html.Div(
                        [
                            html.Button('Executar Fluxo de Potência', id='exe-EE', n_clicks=0, style = {'color':'#FFFFFF' , 'background-color' : '#98151b'},), 
                            
                        ],
                        style = {'text-align':'center' , 'margin-top': '300px'},
                        ),

                        

                    html.Div(
                                [
                                    html.H6("Grafo da Rede Recebida", className="subtitle padded"),
                                    cyto.Cytoscape(
                                        id='cytoscape_FPD',
                                        elements={},
                                        selectedEdgeData=[],
                                        selectedNodeData=[],
                                        layout={'name': 'cose'}, # 'cose','grid'
                                        style={'width': '700px', 'height': '500px', 'border': '1px solid gray', 'border-radius': '5px'},
                                        stylesheet = [
                                            {
                                                'selector': 'node',
                                                'style': {
                                                    'label': 'data(id)',
                                                    'shape': 'polygon',
                                                    'shape-polygon-points': '-0.2, -1, 0.2, -1, 0.2, 1, -0.2, 1',
                                                    'background-color': 'black'
                                                }
                                            },
                                            {
                                                'selector': ':selected',
                                                'style': {
                                                    'background-color': 'SteelBlue',
                                                }
                                            },
                                            {
                                                'selector': 'edge',
                                                'style': {
                                                    # "curve-style": "taxi",
                                                    "taxi-turn": 20,
                                                    "taxi-turn-min-distance": 5
                                                }
        
                                            }
                                        ],
                                        responsive = True
                                    ),
                                ],
                                className="row ",
                            ),        

                    #Gráfico do processo Iterativo
                    html.Div([dcc.Graph(id="graphFPD")]),


                    #Tabela dos Resultados do Fluxo
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H6(
                                        ["Resultado de Medidas das Barras"], className="subtitle padded"
                                    ),
                                    dash_table.DataTable(id = 'FPD_bar',page_action='none',
                                    style_table={'overflow': 'auto', 'width': '500px','align':'center'},style_cell={'textAlign': 'center'}),
                                ],              
                                       
                            ),   
                        ],
                        className="row ",
                    ),

                    html.Div(
                        [
                        html.Button("Download Excel", id="btn_bar_xlsx"),
                        dcc.Download(id="download-FPD_Bar-xlsx"),
                        ],style={'margin-bottom':50}
                    ),

                    #Tabela dos Resultados do Fluxo
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H6(
                                        ["Resultados de Medidas das Linhas"], className="subtitle padded"
                                    ),
                                    dash_table.DataTable(id = 'FPD_lin',page_action='none', #fixed_rows={'headers': True}, 
                                    style_table={'overflow': 'auto', 'width': '600px','align':'center'},style_cell={'textAlign': 'center'}),
                                ],                 
                                       
                            ),   
                        ],
                        className="row ",
                    ),

                    html.Div(
                        [
                        html.Button("Download Excel", id="btn_lin_xlsx"),
                        dcc.Download(id="download-FPD_Lin-xlsx"),
                        ],style={'margin-bottom':50}
                    ),

                    
             ],
             
             className="sub_page",
             ),
             
           
        ],
        className="page",  
    )
