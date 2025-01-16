import sys
import os

# PACKAGE_PARENT = '..'
# SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
# sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

import dash
from networkx.classes.function import is_empty
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State
from dash_table.Format import Format, Group, Scheme, Symbol
from dash_app import app
import multilanguage as ml
from callbacks.comuns import *
import numpy as np
from SisPotFunctions import PSFPD
import dash_core_components as dcc
import traceback




############################################################Power Flow##########################################################################################
@app.callback(
    Output("download-FPD_Lin-xlsx", "data"),
    Input("btn_lin_xlsx", "n_clicks"),
    State("FPD_lin","data"),
    prevent_initial_call=True,
)
def download_res_linha(n_clicks,data):
    df_lin = pd.DataFrame(data)
    if (df_lin.empty):return None
    else: return dcc.send_data_frame(df_lin.to_excel, "Line_Results.xlsx")

@app.callback(
    Output("download-FPD_Bar-xlsx", "data"),
    Input("btn_bar_xlsx", "n_clicks"),
    State("FPD_bar","data"),
    prevent_initial_call=True,
)
def download_res_bar(n_clicks,data):
    df_lin = pd.DataFrame(data)
    if (df_lin.empty):return None
    else: return dcc.send_data_frame(df_lin.to_excel, "Bus_Results.xlsx")

@app.callback(
    Output("download-MedsFPD-xlsx", "data"),
    Input("btn_meds_xlsx", "n_clicks"),
    State("FPD_meds","data"),
    prevent_initial_call=True,
)
def download_res_meds(n_clicks,data):
    df_lin = pd.DataFrame(data)
    if (df_lin.empty):return None
    else: return dcc.send_data_frame(df_lin.to_excel, "Measurement_Results.xlsx")


@app.callback([Output("topology_table_FPD", "data"), 
               Output("topology_table_FPD", "columns"),
               Output('cytoscape_FPD','elements')],
              [Input('topology_FPD', 'contents')],
              [State('topology_FPD', 'filename'),
               State('tables-storage', 'data'),
               State('page-language', 'lang')])
def insert_topology_fpd(contents, filename, tablesData, language):
    if not contents:
        tablesDataHasContent = False
        if tablesData:
            tablesDataHasContent = "TopologyData" in tablesData
        if tablesDataHasContent:
            return tablesData['TopologyData'], tablesData['TopologyColumns'], tablesData['TopologyCytoscapeElements']
        else:
            return  [{" ":" "}], [{"name": " ", "id": " "}], {}
    else:
        dff, rede = parse_contents(contents, filename)
        calls=[{"name": ml.ml(x, language), "id": x} if dff[x].dtypes == object else {"name": ml.ml(x, language), "id": x,'format': Format(precision=4),'type':'numeric'} for x in dff.columns]
        return dff.to_dict('records'),calls,rede

# @app.callback([Output("topology_table_FPD", "data")],
#               [State('topology_FPD', 'filename'),
#                State('tables-storage', 'data'),
#                State('page-language', 'lang')])
# def insert_topology_fpd(contents, filename, tablesData, language):


@app.callback([Output("load_table_FPD", "data"),
               Output("load_table_FPD", "columns")],
              [Input('load_FPD', 'contents')],
              [State('load_FPD', 'filename'),
               State('load-storage-FPD', 'data'),
               State('page-language', 'lang')])
def load_fpd(contents, filename, tablesData, language):
    # if not contents:
    #     return  [{" ":" "}], [{"name": " ", "id": " "}]
    if not contents:
        tablesDataHasContent = False
        if tablesData:
            tablesDataHasContent = "LoadData" in tablesData
        if tablesDataHasContent:
            return tablesData['LoadData'], tablesData['LoadColumns']
        else:
            return  [{" ":" "}], [{"name": " ", "id": " "}]
    else:
        dff, _ = parse_contents(contents, filename,False)
        calls=[{"name": ml.ml(x, language), "id": x} if dff[x].dtypes == object else {"name": ml.ml(x, language), "id": x,'format': Format(precision=4),'type':'numeric'} for x in dff.columns]
        return dff.to_dict('records'),calls

@app.callback([Output("FPD_bar", "data"), Output("FPD_bar", "columns"),Output("FPD_lin", "data"), Output("FPD_lin", "columns"),Output("graphFPD","figure"),Output('confirm-error-fpd','displayed'),Output('confirm-error-fpd', 'message')], 
              [Input("topology_table_FPD", "data"),Input("load_table_FPD", "data"),Input('input_tol_FPD', 'value'),Input('input_iteracoes_FPD', 'value'),Input('input_refbus_FPD', 'value'),Input('input_final_bus_FPD', 'value'),Input('exe-EE', 'n_clicks')],
              State('page-language', 'lang'))
def update_fpd(topology_data,load_data,tol,max_it,ref,final_ref,n_clicks,language):
    
    try:
    
        if len(topology_data) > 2 and len(load_data) > 2:

            network_file = pd.DataFrame(topology_data)
            load_file = pd.DataFrame(load_data).dropna() 

            if len(dash.callback_context.triggered):
                if dash.callback_context.triggered[0]['prop_id'] == 'exe-EE.n_clicks':
                    
                    if ref == final_ref: 
                        return [{" ":" "}], [{"name": " ", "id": " "}],[{" ":" "}],[{"name": " ", "id": " "}],{},True,'A referência deve ser diferente da ultima barra'
                    if not PSFPD.are_endings(ref,topology_data):
                        return [{" ":" "}], [{"name": " ", "id": " "}],[{" ":" "}],[{"name": " ", "id": " "}],{},True,'A barra de Referência deve estar nas extremidades da rede'
                    if not PSFPD.are_endings(final_ref,topology_data):
                        return [{" ":" "}], [{"name": " ", "id": " "}],[{" ":" "}],[{"name": " ", "id": " "}],{},True,'A barra final do ramo principal deve estar nas extremidades da rede'
                    
                    resultados_barra, resultados_fluxo,resultados_linha,resultado_medidas, iter, n_bus,objetivo= PSFPD.run_fpd(network_file,load_file,tol,max_it,ref,final_ref) #roda e armazena o resultado do fluxo de potencia
                    
                    #########################################################Prepara Plotagem#######################################################
                    ativo = [x[0] for x in objetivo]
                    reativo = [x[1] for x in objetivo]
                    fig = px.line(x=list(range(1,len(ativo)+1)),y=ativo,title=ml.ml('Desvio Máximo por Iteração', language),labels=dict(x=ml.ml('Iterações', language),y=ml.ml('Desvio', language),name=ml.ml('Desvio Ativo', language))) 
                    
                    
                    fig.add_trace(go.Scatter(x=list(range(1, len(ativo)+1))  , y=ativo  , mode='lines', name=ml.ml('Desvio Ativo'  ,language)))
                    fig.add_trace(go.Scatter(x=list(range(1, len(reativo)+1)), y=reativo, mode='lines', name=ml.ml('Desvio Reativo',language)))
                    fig.add_hline(y=0.0001, line_width=3, line_dash="dash", line_color="green")
                    fig.update_yaxes(exponentformat="power")
                    fig.update_xaxes(dtick=1)
                    
                    # #####################################################################################################
                    
                    cols=[{"name": ml.ml(x, language), "id": x} if resultados_barra[x].dtypes == object else {"name": ml.ml(x, language), "id": x,'format': Format(precision=4),'type':'numeric'} for x in resultados_barra.columns]
                    cols2=[{"name": ml.ml(x, language), "id": x} if resultados_fluxo[x].dtypes == object else {"name": ml.ml(x, language), "id": x,'format': Format(precision=4),'type':'numeric'} for x in resultados_fluxo.columns]
                    # calls3=[{"name": ml.ml(x, language), "id": x} if resultados_linha[x].dtypes == object else {"name": ml.ml(x, language), "id": x,'format': Format(precision=4),'type':'numeric'} for x in resultados_linha.columns]
                    return resultados_barra.to_dict('records'),cols,resultados_fluxo.to_dict('records'),cols2,fig,True,'Fluxo calculado com sucesso!'

        return [{" ":" "}], [{"name": " ", "id": " "}],[{" ":" "}],[{"name": " ", "id": " "}],{},False,''
    
    except IndexError as e:
        print(f'Erro de dimensão!: {e}')
        return [{" ":" "}], [{"name": " ", "id": " "}],[{" ":" "}],[{"name": " ", "id": " "}],{},True,'ERRO na dimensão dos arquivos de entrada'
    except Exception as e:
        print(f'Erro inesperado!: {e}')
        traceback.print_exc()
        return [{" ":" "}], [{"name": " ", "id": " "}],[{" ":" "}],[{"name": " ", "id": " "}],{},True,'ERRO INESPERADO! Entre em contato com o suporte!'
        

###### global system data updater routines ######
@app.callback(
    Output('load-storage-FPD', 'data'),
    Input("load_table_FPD", "data"),
    Input("load_table_FPD", "columns"), 
    State('load-storage-FPD', 'data'),
)
def saveTablesFP(topologyData, topologyColumns, tablesData):
    # if ts is None:
    #     raise PreventUpdate
    if not tablesData:
        tablesData = dict()
    tablesData['LoadData'] = topologyData
    tablesData['LoadColumns'] = topologyColumns
    return tablesData

