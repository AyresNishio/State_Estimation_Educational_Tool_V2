import numpy as np
import pandas as pd
import os
from scipy.stats import chi2

from SisPotFunctions.PF_functions import *


def state_estimation( line_data, meas_data, tol, metodo,max_it):
    tol = float(tol)
    max_it = int(max_it)
    Y_bus, network_values = build_Y_bus(line_data)
    n_bus = len(Y_bus)
    for i in range(n_bus):
        for j in range(n_bus):
            
            n = Y_bus[i,j]
            im='({c.real:.4f} + {c.imag:.4f}i)'.format(c=n)
            print(im,end = ' ')
        print()
    meas, R, order_meas = read_measures(meas_data)

    is_observable, criticality_data = linear_SE(order_meas, Y_bus)
    
    if (is_observable):
        #print('observable')
        if  (metodo == 'padrao'):
            J_list,J_critical,State, Filtered_meas,Saidas = standard_SE(Y_bus, order_meas, meas, network_values,R,tol,max_it)
        elif(metodo == 'algoritmo'):
            J_list,J_critical,State, Filtered_meas,Saidas = decoupled_SE_algorithm(Y_bus, order_meas, meas, network_values,R,tol,max_it)
        elif(metodo == 'modelo'):
            J_list,J_critical,State, Filtered_meas,Saidas = decoupled_SE_model(Y_bus, order_meas, meas, network_values,R,tol,max_it)
        return criticality_data, J_list, J_critical, State, Filtered_meas, is_observable, Saidas    
    else:
        return 0, 0, 0, 0, 0, is_observable, 0

def build_Y_bus( line_data):
    
    #path_file = os.path.join(filepath, filename)
    #line_data = pd.read_excel(path_file)
    #dict_names = {'De': 'De','Para':'Para','R':'R','X':'X','B':'C','Tap':'Tap'}
    #line_data.rename(columns=dict_names,inplace=True)
    
    total_lines = set(line_data['De'].to_list() + line_data['Para'].to_list())
    n = len(total_lines)
    Y_bus = np.zeros((n, n), dtype=complex)
    connective_values = {}
    for lines in line_data.iterrows():
        data_lines = lines[1]
        bus_from = int(data_lines['De'] - 1)
        bus_to = int(data_lines['Para'] - 1)
        coord = (bus_from, bus_to)
        bshunt = float(data_lines['B'])
        diag_value, out_diag_value = get_line_params(data_lines['R'], data_lines['X'], bshunt, data_lines['Tap'])
        _, line_admitance = get_line_params(data_lines['R'], data_lines['X'], 0, 0)
        line_admitance_shunt = bshunt*0.5j if bshunt > 0 else 0
        connect = {('g', *coord): line_admitance.real, 
                   ('b', *coord): line_admitance.imag,
                   ('gs', *coord): line_admitance_shunt.real if type(line_admitance_shunt) == complex else 0, 
                   ('bs', *coord): line_admitance_shunt.imag if type(line_admitance_shunt) == complex else 0}
        if len(connective_values) == 0:
            connective_values = connect.copy()
        else:
            connective_values.update(connect)
        Y_bus[bus_from, bus_from] += diag_value
        Y_bus[bus_to, bus_to] += diag_value if data_lines['Tap'] == 0 else diag_value*(data_lines['Tap']**2)
        Y_bus[bus_from, bus_to] -= out_diag_value
        Y_bus[bus_to, bus_from] -= out_diag_value
        if bus_to ==1 or bus_from ==1: print(diag_value)
    return Y_bus, connective_values    

def get_line_params(r: float, x: float, c_total: float, a: float) -> complex:
    diag_value = 1/(r + x*1j) + c_total*0.5j if a == 0 else (1/(r + x*1j) + c_total*0.5j)/(a**2)
    out_diag_value = 1/(r + x*1j) if a == 0 else (1/(r + x*1j))/(a)
    return diag_value, out_diag_value

def read_measures(measured_values):
       
    measures = {}
    R = np.identity(measured_values.shape[0])
    order_meansured = []
    for i, line in enumerate(measured_values.iterrows()):
        meansured_line = line[1]
        meansured_type = meansured_line['Tipo']
        bus_from_to = (meansured_line['De'], meansured_line['Para']) if meansured_line['Para'] != '-' else meansured_line['De']
        if meansured_type == "P" or meansured_type == "Q" or meansured_type == "I":
            meansured_type_specific = f'{meansured_type}_flow' if type(bus_from_to) == tuple else f'{meansured_type}_inj'
        else:
            meansured_type_specific = meansured_type
        order_meansured.append({meansured_type_specific: bus_from_to})
        R[i,i] = meansured_line['Desvio Padrão']**2
        med_caracteristic = (meansured_line['Valor'], meansured_line['Desvio Padrão'])
        data_med = {bus_from_to: med_caracteristic}
        if meansured_type not in list(measures.keys()):
            measures[meansured_type] = data_med
        else:
            measures[meansured_type].update(data_med)
    return measures, R, order_meansured

#--Observabilidade
def linear_SE(meas, Ybus):
    med_data = pd.DataFrame(columns = ['Tipo', 'Localização'])

    #Linear State Estimation
    line_amount = np.sum([1 if list(x.keys())[0][0] == 'P' else 0 for x in meas ])
    col_amount = Ybus.shape[0] - 1
    H_linear = np.zeros((line_amount, col_amount))
    l = 0
    for meansured in meas:
        meansured_type = list(meansured.keys())[0]
        meansured_loc = list(meansured.values())[0]
        if meansured_type[0] != 'P':
            pass
        else:
            #med_data.append([meansured_type, meansured_loc])
            med_data.loc[len(med_data)] = [meansured_type, meansured_loc]
            if meansured_type.split('_')[-1] == 'flow':
                if meansured_loc[0] - 2 >= 0:
                    H_linear[l][meansured_loc[0] - 2] = 1
                if meansured_loc[1] - 2 >= 0:
                    H_linear[l][meansured_loc[1] - 2] = -1
            else:
                aux = Ybus[meansured_loc - 1][1:]
                conected = aux == 0
                H_linear[l][~conected] = -1
                if meansured_loc - 2 >= 0:
                    H_linear[l][meansured_loc - 2] = np.sum(conected)
            l += 1
    R = np.identity(line_amount)
    G = H_linear.T@H_linear
    try:
        detG = np.linalg.det(G)
        G_1 = np.linalg.inv(G)
        E = R - (H_linear@(G_1))@H_linear.T
    except:
        detG = 0
        G_1 = np.zeros(G.shape)
        E = np.zeros(R.shape)

    #Observability Test
    is_observable=False if detG <= 1e-11 else True
    criticality_data = []
    if is_observable : 
        criticality_data = identify_critical(E,med_data)
    
    return is_observable, criticality_data

#--Estimação de Estado
class dados_iter():
    def __init__(self,hx,Jacob,G_matrix,delta_x,state_array) -> None:
        self.h_x = hx.copy()
        self.Jacob = Jacob.copy()
        self.G_matrix = G_matrix.copy()
        self.delta_x = delta_x.copy()
        self.state_array = state_array.copy()
class dados_decoupled_iter():        
    def __init__(self,hx1,hx2,Jacob_T,Gain_T,Jacob_V,Gain_V,delta_x,state_array) -> None:
        self.h_x1 = hx1.copy()
        self.h_x2 = hx2.copy()
        self.Jacob_T = Jacob_T.copy()
        self.Gain_T = Gain_T.copy()
        self.Jacob_V = Jacob_V.copy()
        self.Gain_V = Gain_V.copy()
        self.delta_x = delta_x.copy()
        self.state_array = state_array.copy()

#--Metodos        
def standard_SE(Y_bus,order_meas, meas, network_values,R,tol,max_it):
    
    error = tol+1
    

    n_bus = Y_bus.shape[0]
    state = [0]*(n_bus-1)
    state.extend([1]*n_bus)
    state_array = np.array(state, dtype=float)

    nits = 0
    J_list = []
    Saidas = []

    while error > tol and nits<max_it:
        nits += 1

        h_x, estimated_measures = filter_measurements(order_meas, meas, Y_bus, network_values, state_array)
        H = build_Jacobian(order_meas, Y_bus, network_values, state_array)
        
        G = build_Gain(H, R)
        H_t = np.transpose(H)
        R_inv = np.linalg.inv(R)

        h_x = np.where(np.abs(h_x) < 1e-6, 0, h_x)

        J = np.sum(np.dot(h_x**2, R_inv))
        J_list.append(J)

        h_r = H_t@R_inv
        t = h_r@h_x

        delta_x = np.linalg.solve(G, t)
        state_array += delta_x
        error = np.max(np.abs(delta_x))

        Saidas.append(dados_iter(estimated_measures,H,G,delta_x,state_array))


    #Data_frame dos estado
    State_dataframe = pd.DataFrame(state_array[n_bus-1:], columns = ['Mag. de Tensão'], index = list(range(1,n_bus+1)))
    State_dataframe['Ang.(°)'] = np.degrees([0]+state_array[:n_bus-1].tolist())


    #Data_frame das medidas
    meas_dataframe = pd.DataFrame([med.keys() for med in order_meas], columns = ['Tipos'])
    meas_dataframe['Localização'] = [list(loc.values())[0] for loc in order_meas]
    meas_dataframe['Valor Medido'] = [meas[list(x.keys())[0].split('_')[0]][list(x.values())[0]][0] for x in order_meas]
    meas_dataframe['Valor Estimado'] = estimated_measures
    meas_dataframe['Desvio'] = np.where(np.abs(h_x) < 1e-6, 0, h_x)
    Sigma = np.linalg.inv(G)
    Erii = np.abs(R - H@Sigma@H_t)
    Erii = np.sqrt(np.diag(Erii))
    meas_dataframe['Res. Normalizado'] = np.where(Erii > 10e-10, np.abs(meas_dataframe['Desvio']/Erii), np.inf)
    meas_dataframe['Tipos'] = meas_dataframe['Tipos'].map({"P_flow": "Fluxo Pot. Ativ.", "Q_flow": "Fluxo Pot. Reativ.", "P_inj": "Injeção Pot. Ativ.", "Q_inj": "Injeção Pot. Reativ.",'V':'Módulo da Tensão'})

    #J_critico
    ddof = meas_dataframe.shape[0] - 2*(State_dataframe.shape[0]) + 1
    J_critical = chi2.ppf(1-0.05, ddof)

    texto_saida  = print_steps(Saidas,State_dataframe,meas_dataframe)
    return J_list,J_critical,State_dataframe, meas_dataframe,texto_saida

def decoupled_SE_algorithm(Y_bus,order_meas, meas, network_values,R,tol,max_it):
       
    error = tol+1
    
    n_bus = Y_bus.shape[0]
    state = [0]*(n_bus-1)
    state.extend([1]*n_bus)
    state_array = np.array(state, dtype=float)

    medidas_ativas = []
    medidas_reativas = []
    for meas2 in order_meas:
        if 'P_flow' in meas2 or 'P_inj' in meas2: medidas_ativas.append(meas2)
        if 'Q_flow' in meas2 or 'Q_inj'in meas2 or 'V' in meas2: medidas_reativas.append(meas2)
    order_meas = []
    order_meas.extend(medidas_ativas)
    order_meas.extend(medidas_reativas)
    range_reativas = (len(medidas_ativas),len(order_meas)-1)
    
    nits = 0
    J_list = []
    Saidas = []
    
    nvs=int(np.floor((len(state))/2))

    while error > tol and nits<max_it:
        nits += 1

        h_x, estimated_measures = filter_measurements(order_meas, meas, Y_bus, network_values, state_array)
        H = build_Jacobian(order_meas, Y_bus, network_values, state_array)
        for i in  range(0,range_reativas[1]+1):
            for j in range(0,len(state)):
                if i<len(medidas_ativas):
                    if j>=nvs:
                        H[i,j] = 0
                if j<nvs:
                    if i>=len(medidas_ativas):
                        H[i,j] = 0
        G = build_Gain(H, R)
                
        Gt = G[0:nvs,0:nvs]
        
        Ht = H[0:len(medidas_ativas),0:nvs]
        
        Rt = R[0:len(medidas_ativas),0:len(medidas_ativas)]
                
        delta_z, estimated_measures = filter_measurements(order_meas, meas, Y_bus, network_values, state_array)
        e_values1 = estimated_measures[:len(medidas_ativas)]
        delta_z = np.where(np.abs(delta_z) < 1e-6, 0, delta_z)
        
        zt = delta_z[0:len(medidas_ativas)]

        Ht_t = np.transpose(Ht)
        Rt_inv = np.linalg.inv(Rt)
        
        zt_r = Ht_t@Rt_inv
        tt = zt_r@zt

        delta_xt = np.linalg.solve(Gt, tt)
        state_array[0:nvs] += delta_xt
        error_t = np.max(np.abs(delta_xt))
        
        H = build_Jacobian(order_meas, Y_bus, network_values, state_array)
        for i in  range(0,range_reativas[1]+1):
            for j in range(0,len(state)):
                if i<len(medidas_ativas):
                    if j>=nvs:
                        H[i,j] = 0
                if j<nvs:
                    if i>=len(medidas_ativas):
                        H[i,j] = 0               
        G = build_Gain(H, R)
        
        Gv = G[nvs:,nvs:]
        
        Hv = H[len(medidas_ativas):,nvs:]
        
        Rv = R[len(medidas_ativas): ,len(medidas_ativas):]
        
        

        delta_z, estimated_measures = filter_measurements(order_meas, meas, Y_bus, network_values, state_array)
        e_values2 = estimated_measures[len(medidas_ativas):]
        delta_z = np.where(np.abs(delta_z) < 1e-6, 0, delta_z)
        
        zv = delta_z[len(medidas_ativas):]

        Hv_t = np.transpose(Hv)
        Rv_inv = np.linalg.inv(Rv)
        zv_r = Hv_t@Rv_inv
        tv = zv_r@zv

        delta_xv = np.linalg.solve(Gv, tv)
        state_array[nvs:] += delta_xv
        error_v = np.max(np.abs(delta_xv))
        
        error = max(error_t,error_v)
        delta_x = np.append(delta_xt,delta_xv)
        
        Saidas.append(dados_decoupled_iter(e_values1,e_values2,Ht,Gt,Hv,Gv,delta_x,state_array))
        H_t = np.transpose(H)
        R_inv = np.linalg.inv(R)
        J = np.sum(np.dot(delta_z**2, R_inv))
        J_list.append(J)
        


    #Data_frame dos estado
    State_dataframe = pd.DataFrame(state_array[n_bus-1:], columns = ['Mag. de Tensão'], index = list(range(1,n_bus+1)))
    State_dataframe['Ang.(°)'] = np.degrees([0]+state_array[:n_bus-1].tolist())

    #Data_frame das medidas
    meas_dataframe = pd.DataFrame([med.keys() for med in order_meas], columns = ['Tipos'])
    meas_dataframe['Localização'] = [list(loc.values())[0] for loc in order_meas]
    meas_dataframe['Valor Medido'] = [meas[list(x.keys())[0].split('_')[0]][list(x.values())[0]][0] for x in order_meas]
    meas_dataframe['Valor Estimado'] = estimated_measures
    meas_dataframe['Desvio'] = np.where(np.abs(h_x) < 1e-6, 0, h_x)
    Sigma = np.linalg.inv(G)
    Erii = np.abs(R - H@Sigma@H_t)
    Erii = np.sqrt(np.diag(Erii))
    meas_dataframe['Res. Normalizado'] = np.where(Erii > 10e-10, np.abs(meas_dataframe['Desvio']/Erii), np.inf)
    meas_dataframe['Tipos'] = meas_dataframe['Tipos'].map({"P_flow": "Fluxo Pot. Ativ.", "Q_flow": "Fluxo Pot. Reativ.", "P_inj": "Injeção Pot. Ativ.", "Q_inj": "Injeção Pot. Reativ.",'V':'Módulo da Tensão'})

    #J_critico
    ddof = meas_dataframe.shape[0] - 2*(State_dataframe.shape[0]) + 1
    J_critical = chi2.ppf(1-0.05, ddof)
    
    texto_saida  = print_steps_decoupled(Saidas,State_dataframe,meas_dataframe)
    return J_list,J_critical,State_dataframe, meas_dataframe,texto_saida


def decoupled_SE_model(Y_bus,order_meas, meas, network_values,R,tol,max_it):
       
    error = tol+1
    
    n_bus = Y_bus.shape[0]
    state = [0]*(n_bus-1)
    state.extend([1]*n_bus)
    state_array = np.array(state, dtype=float)

    medidas_ativas = []
    medidas_reativas = []
    for meas2 in order_meas:
        if 'P_flow' in meas2 or 'P_inj' in meas2: medidas_ativas.append(meas2)
        if 'Q_flow' in meas2 or 'Q_inj'in meas2 or 'V' in meas2: medidas_reativas.append(meas2)
    order_meas = []
    order_meas.extend(medidas_ativas)
    order_meas.extend(medidas_reativas)
    
    nits = 0
    J_list = []
    Saidas = []

    while error > tol and nits<max_it:
        nits += 1

        
        H = build_Jacobian(order_meas, Y_bus, network_values, state_array)
                        
        G = build_Gain(H, R)
        
        nvs=int(np.floor((len(state))/2))
        
        for i in range(len(G)):
            for j in range(len(G[0])):
                if i<nvs:
                    if j>=nvs:
                        G[i,j] = 0
                if j<nvs:
                    if i>=nvs:
                        G[i,j] = 0
        
        Gt = G[0:nvs,0:nvs]
        
        Ht = H[0:len(medidas_ativas),0:nvs]
        
        Rt = R[0:len(medidas_ativas),0:len(medidas_ativas)]
                
        delta_z, estimated_measures = filter_measurements(order_meas, meas, Y_bus, network_values, state_array)
        e_values1 = estimated_measures[:len(medidas_ativas)]
        delta_z = np.where(np.abs(delta_z) < 1e-6, 0, delta_z)
        
        zt = delta_z[0:len(medidas_ativas)]

        Ht_t = np.transpose(Ht)
        Rt_inv = np.linalg.inv(Rt)
        
        zt_r = Ht_t@Rt_inv
        tt = zt_r@zt

        delta_xt = np.linalg.solve(Gt, tt)
        state_array[0:nvs] += delta_xt
        error_t = np.max(np.abs(delta_xt))
        
        H = build_Jacobian(order_meas, Y_bus, network_values, state_array)
                        
        G = build_Gain(H, R)
        for i in range(len(G)):
            for j in range(len(G[0])):
                if i<nvs:
                    if j>=nvs:
                        G[i,j] = 0
                if j<nvs:
                    if i>=nvs:
                        G[i,j] = 0
        
        Gv = G[nvs:,nvs:]
        
        Hv = H[len(medidas_ativas):,nvs:]
        
        Rv = R[len(medidas_ativas): ,len(medidas_ativas):]
        
        

        delta_z, estimated_measures = filter_measurements(order_meas, meas, Y_bus, network_values, state_array)
        e_values2 = estimated_measures[len(medidas_ativas):]
        delta_z = np.where(np.abs(delta_z) < 1e-6, 0, delta_z)
        
        zv = delta_z[len(medidas_ativas):]

        Hv_t = np.transpose(Hv)
        Rv_inv = np.linalg.inv(Rv)
        zv_r = Hv_t@Rv_inv
        tv = zv_r@zv

        delta_xv = np.linalg.solve(Gv, tv)
        state_array[nvs:] += delta_xv
        error_v = np.max(np.abs(delta_xv))
        
        error = max(error_t,error_v)
        delta_x = np.append(delta_xt,delta_xv)
        
        Saidas.append(dados_decoupled_iter(e_values1,e_values2,Ht,Gt,Hv,Gv,delta_x,state_array))
        H_t = np.transpose(H)
        R_inv = np.linalg.inv(R)
        J = np.sum(np.dot(delta_z**2, R_inv))
        J_list.append(J)

    

    #Data_frame dos estado
    State_dataframe = pd.DataFrame(state_array[n_bus-1:], columns = ['Mag. de Tensão'], index = list(range(1,n_bus+1)))
    State_dataframe['Ang.(°)'] = np.degrees([0]+state_array[:n_bus-1].tolist())


    #Data_frame das medidas
    meas_dataframe = pd.DataFrame([med.keys() for med in order_meas], columns = ['Tipos'])
    meas_dataframe['Localização'] = [list(loc.values())[0] for loc in order_meas]
    meas_dataframe['Valor Medido'] = [meas[list(x.keys())[0].split('_')[0]][list(x.values())[0]][0] for x in order_meas]
    meas_dataframe['Valor Estimado'] = estimated_measures
    meas_dataframe['Desvio'] = np.where(np.abs(delta_z) < 1e-6, 0, delta_z)
    Sigma = np.linalg.inv(G)
    Erii = np.abs(R - H@Sigma@H_t)
    Erii = np.sqrt(np.diag(Erii))
    meas_dataframe['Res. Normalizado'] = np.where(Erii > 10e-10, np.abs(meas_dataframe['Desvio']/Erii), np.inf)
    meas_dataframe['Tipos'] = meas_dataframe['Tipos'].map({"P_flow": "Fluxo Pot. Ativ.", "Q_flow": "Fluxo Pot. Reativ.", "P_inj": "Injeção Pot. Ativ.", "Q_inj": "Injeção Pot. Reativ.",'V':'Módulo da Tensão'})

    #J_critico
    ddof = meas_dataframe.shape[0] - 2*(State_dataframe.shape[0]) + 1
    J_critical = chi2.ppf(1-0.05, ddof)

    texto_saida  = print_steps_decoupled(Saidas,State_dataframe,meas_dataframe)
    return J_list,J_critical,State_dataframe, meas_dataframe,texto_saida

#--Etapas da EE
def filter_measurements(order_meansured, meansured_network, Y_bus, network_values, state_vetctor):
    n_bus = Y_bus.shape[0]
    if state_vetctor is None:
        states = [0]*(n_bus-1)
        states.extend([1]*n_bus)
        state_vetctor = np.array(states)

    n_lines = len(order_meansured)
    z_array = np.zeros(shape = (n_lines, ))
    estimated_measurements = np.zeros(shape = (n_lines, ))
    func_dict = {'P_inj': p_Inj,'P_flow': p_Flow, 'Q_inj': q_Inj, 'Q_flow': q_Flow, 'V': voltage}
    
    for i, meansured in enumerate(order_meansured):
        meansured_type = list(meansured.keys())[0]
        bus = meansured.get(meansured_type)
        line_h = func_dict[meansured_type](Y_bus, network_values, state_vetctor, bus, n_bus)
        estimated_measurements[i] = line_h
        z_array[i] = meansured_network[meansured_type[0]][bus][0] - line_h

    return z_array, estimated_measurements

def build_Jacobian(order_meansured: list, y_bar: np.ndarray, network_values: dict, state_vetctor: np.ndarray = None):
    n_bus = y_bar.shape[0]
    if state_vetctor is None:
        states = [0]*(n_bus-1)
        states.extend([1]*n_bus)
        state_vetctor = np.array(states)

    line_amount = len(order_meansured)
    col_amount  = len(state_vetctor)
    h_matrix = np.zeros(shape = (line_amount, col_amount))
    func_dict = {'P_inj': derivada_p_Inj,'P_flow': derivada_p_Flow, 'Q_inj': derivada_q_Inj, 'Q_flow': derivada_q_Flow, 'V': derivada_voltage}
    
    for i, meansured in enumerate(order_meansured):
        meansured_type = list(meansured.keys())[0]
        bus = meansured.get(meansured_type)
        line_h = func_dict[meansured_type](y_bar, network_values, state_vetctor, bus, n_bus)
        h_matrix[i,:] = line_h
    return h_matrix

def build_Gain(H: np.ndarray, R: np.ndarray) -> np.ndarray:
    R_inv = np.linalg.inv(R)
    H_T = np.transpose(H.copy())
    h_r_inv = H_T@R_inv
    G = h_r_inv@H
    return G

#-- Analise de criticalidade
def identify_critical(E: np.ndarray, med_data: pd.core.frame.DataFrame):
    med_data['Criticalidades'] = np.nan
    amount_meansured = E.shape[0]
    residual = np.sum(E, axis = 1)
    standart_dev = np.sqrt(np.diag(E))
    aux = standart_dev[:,None]*standart_dev
    aux = np.where(aux == 0, 0.00000001, aux)
    gamma = np.divide(np.abs(E), aux)
    normalized_residual = np.abs(np.where(standart_dev == 0, 0, residual/standart_dev))
    normalized_residual__ = np.where(normalized_residual == 0, 0.000001, normalized_residual)
    rho = normalized_residual__[:,None]/normalized_residual__
    critical_meansured = (standart_dev <= 1e-6)*(residual <= 1e-6)
    number_cmeans = list(np.arange(amount_meansured)[critical_meansured])
    med_data.loc[med_data.index.isin(number_cmeans), 'Criticalidades'] = 'Medida Crítica'
    critical_sets = []
    csets = []
    amount_csets = 0
    for i in range(0, amount_meansured):
        for j in range(i+1, amount_meansured):
            if (rho[i][j] >= 0.98) and (gamma[i][j] >= 0.98) and (i not in number_cmeans) and (j not in number_cmeans):
                csets.append(i)
                csets.append(j)
        if any([set(csets).issubset(x) for x in critical_sets]):
            csets = []
        if len(csets) > 0:
            amount_csets += 1
            critical_sets += [set(csets)]
            med_data.loc[med_data.index.isin(csets), 'Criticalidades'] = f'Conj.Crítico_{amount_csets}'
            csets = []
    return med_data


#--Saidas .txt
def print_steps(saidas_it, state_dataframe, meas_dataframe):
    output_text = ''

    it = 0
    for saida in saidas_it:

        h_x         = saida.h_x 
        Jacob       = saida.Jacob 
        G_matrix    = saida.G_matrix 
        delta_x     = saida.delta_x 
        state_array = saida.state_array  

        output_text += f'\n##### ITERATION {it+1} #####\n'
        output_text += '\n### Jacobian Matrix H ###\n'
        for lin in Jacob:
            for x in lin:
                output_text += format(x,' .4e') + ' '
            output_text+='\n'

        output_text += '\n### Gain Matrix G ###\n'
        for lin in G_matrix:
            for x in lin:
                output_text += format(x,' .4e') + ' '
            output_text+='\n'
        
        output_text += '\n### Estimated Values h(x) ###\n'
        for x in h_x:
            output_text += format(x,' .4f')+'\n'
        output_text += '\n' 

        output_text += '\n### Delta x ###\n'
        for x in delta_x:
            output_text += format(x,' .10f')+'\n'
        output_text += '\n' 

        output_text += '\n### State (x) ###\n'
        for x in state_array:
            output_text += format(x,' .4f')+'\n'
        output_text += '\n'

        it+=1

    output_text += f'\n##### Final Results #####\n'
    output_text += '\n### Final Estimated State ###\n'
    output_text += state_dataframe.to_string() +'\n'

    #output_text += '\n### Valores Medidos Finais###\n'
    output_text += '\n### Final Estimated Values ###\n'
    output_text += meas_dataframe.to_string()
    return output_text

def print_steps_decoupled(saidas_it, state_dataframe, meas_dataframe):
    output_text = ''

    it = 0
    for saida in saidas_it:

        h_x1         = saida.h_x1 
        h_x2         = saida.h_x2
        Jacob_t       = saida.Jacob_T
        Jacob_v       = saida.Jacob_V
        Gain_t    = saida.Gain_T
        Gain_v    = saida.Gain_V
        delta_x     = saida.delta_x 
        state_array = saida.state_array  

        output_text += f'\n##### ITERATION {it+1} #####\n'
        output_text += '\n### Jacobian Matrix (Htheta) : Ptheta ###\n'
        for lin in Jacob_t:
            for x in lin:
                output_text += format(x,' .4e') + ' '
            output_text+='\n'
        output_text += '\n### Gain Matrix (Gtheta) : Ptheta  ###\n'
        for lin in Gain_t:
            for x in lin:
                output_text += format(x,' .4e') + ' '
            output_text+='\n'  
        output_text += '\n### Estimated Values h(x) : Ptheta ###\n'
        for x in h_x1:
            output_text += format(x,' .4f')+'\n'
        output_text += '\n'   
            
        output_text += '\n### Jacobian Matrix (HV) : QV ###\n'
        for lin in Jacob_v:
            for x in lin:
                output_text += format(x,' .4e') + ' '
            output_text+='\n'
        output_text += '\n### Gain Matrix (Gv) : QV ###\n'
        for lin in Gain_v:
            for x in lin:
                output_text += format(x,' .4e') + ' '
            output_text+='\n'
        output_text += '\n### Estimated Values h(x) : QV ###\n'
        for x in h_x2:
            output_text += format(x,' .4f')+'\n'
        output_text += '\n' 

        output_text += '\n### Delta x ###\n'
        for x in delta_x:
            output_text += format(x,' .10f')+'\n'
        output_text += '\n' 

        output_text += '\n### State (x) ###\n'
        for x in state_array:
            output_text += format(x,' .4f')+'\n'
        output_text += '\n'

        it+=1

    output_text += f'\n##### FINAL RESULT #####\n'
    output_text += '\n### Final Estimated State ###\n'
    output_text += state_dataframe.to_string() +'\n'

    output_text += '\n### Final Estimated Values ###\n'
    output_text += meas_dataframe.to_string()
    return output_text

if __name__ == "__main__":

    metodo = 'padrao'

    #saidas_it, state_dataframe, meas_dataframe = 
    criticality_data, J_list,J_critical, State_dataframe, data_SE,is_obs=state_estimation('TopologiaCap2.xlsx', 'MedidasCap2.xlsx',0.0001,metodo)
    arquivo = open('Saidas.txt', 'w', encoding='utf-8')
    
    # # Escreve a string no arquivo
    # if(metodo == 'padrao'):
    #     arquivo.write(print_steps(saidas_it,state_dataframe,meas_dataframe))    
    # else:
    #     arquivo.write(print_steps_decoupled(saidas_it,state_dataframe,meas_dataframe))    
    # arquivo.close()
    
    # Saída da Ferramenta _, criticality_data, J_list,J_critical, State_dataframe, data_SE,is_observable
    # print(data_SE)
    # #Pode salvar dessa forma tanto o State_dataframe quanto o data_SE
    # filename = 'Estado_Estimado.csv'
    # State_dataframe.to_csv(filename, sep = ';', decimal = ',')