import pandas as pd
import networkx as nx
from collections import deque
import numpy as np

def are_endings(bus, topology):
    bus=int(bus)
    
    des = [ x['De'] for x in topology]
    paras = [ x['Para'] for x in topology]
    
    grau_ref = des.count(bus) + paras.count(bus) 
    
    if grau_ref>1 : return False
        
    return True
    

def obter_nivel_ramificacao(adj, origem):
    q = deque()
    
    visited = [False] * len(adj)

    visited[origem] = True
    q.append(origem)
    vertices_nivel_lateral = []
    while q:
      
        curr = q.popleft()

        total_adjascencias = adj[curr]
        if (len(total_adjascencias) > 2): 
            vertices_nivel_lateral.append(curr)     

        for x in adj[curr]:
            if not visited[x]:
                visited[x] = True
                q.append(x)

    return vertices_nivel_lateral

def adiciona_arestas(adj, de, para):
    adj[de-1].append(para-1)
    adj[para-1].append(de-1)


def ler_dados(topology_data ,load_data, ref, final_bus):
    param = np.array(topology_data) 
    cargas = np.array(load_data)
    nbus = cargas.shape[0]

    # Processamento das cargas
    carga = {'ativ': [], 'reat': []}
    for i in range(cargas.shape[0]):
        carga['ativ'].append(cargas[i, 1]-cargas[i,3])
        carga['reat'].append(cargas[i, 2]-cargas[i,4])

    matriz_adj = np.zeros((nbus, nbus))
    adj = [[] for _ in range(nbus)]
    G = nx.Graph()
    # Processamento dos parâmetros da rede
    # line = {'de': [], 'para': [], 'r': {}, 'x': {}, 'bsh': {}}
    line = { 'r': {}, 'x': {}, 'bsh': {}}
    for i in range(param.shape[0]):
        prev_bus = int(param[i, 0]) - 1  # Ajustando índice para base 0 (Python)
        next_bus = int(param[i, 1]) - 1

        matriz_adj[prev_bus][next_bus] = 1
        matriz_adj[next_bus][prev_bus] = 1
        G.add_edge(prev_bus, next_bus)
        adiciona_arestas(adj, prev_bus+1, next_bus+1)

        # line['de'].append(prev_bus)
        # line['para'].append(next_bus)
        if (prev_bus, next_bus) not in line['r']:
            line['r'][(prev_bus, next_bus)]   = param[i, 2]
            line['x'][(prev_bus, next_bus)]   = param[i, 3]
            line['bsh'][(prev_bus, next_bus)] = param[i, 4]

    #Processamento Indexação Ramos
    vertices_com_ramificacao = obter_nivel_ramificacao(adj, int(ref)-1) 
    caminho_principal = set(nx.shortest_path(G,int(ref)-1, int(final_bus)-1))
    dicionario_niveis = {}
    for v in vertices_com_ramificacao:
        if v in caminho_principal:
            dicionario_niveis[v+1] = 1
        elif v not in caminho_principal:
            dicionario_niveis[v+1] = 2
    dados = []        
    for bus in range(nbus):
        
        dados.append([])
        dados[bus].append(bus+1)
        
        if bus+1 == int(ref):
            prev_bus = 0
            next_bus = adj[bus][0]+1
        elif G.degree(bus) == 1:
            prev_bus = adj[bus][0]+1
            next_bus  = 0
        elif G.degree(bus) >1:
            prev_bus = adj[bus][0]+1
            next_bus = adj[bus][1]+1
        
        dados[bus].append(prev_bus)
        dados[bus].append(next_bus)
        
        if G.degree(bus) >2:
            dados[bus].append(dicionario_niveis[bus+1])
            barra_ramificada =  adj[bus][-1]+1
            dados[bus].append(barra_ramificada)
            
    return dados, carga, line, nbus, param

def indexacao(dados, nbus):
    # Inicializações
    index = np.zeros((nbus, 4), dtype=int)  # Matriz index com 4 colunas
    maux = np.zeros(nbus, dtype=int)        # Vetor auxiliar
    bf = []                                 # Vetor Breadth-First (BF)
    rbf = []                                # Vetor Reverse Breadth-First (RBF)
    nbf = 0                                 # Número de laterais

    # Configuração inicial
    for i in range(nbus):
        index[i, 0] = dados[i][0]  # Indexação das barras (número da barra na posição 1)

    # Indexação de barras
    for i in range(nbus):
        if i == 0:  # Barra 1 - subestação
            index[i, 1] = 1  # Nível da lateral
            index[i, 2] = 1  # Profundidade da lateral
            index[i, 3] = 0  # Profundidade da barra
            lant = 1         # Nível anterior
            maux[lant - 1] = 1  # Ajuste para índice Python
            nant = 0
        else:  # Demais barras
            if index[i, 1] == 0:  # Nível da barra ainda não determinado
                index[i, 1] = lant  # Atribui nível da lateral ou sublateral
                index[i, 2] = index[int(dados[i][1] - 1), 2]  # Profundidade da lateral (ajuste índice)
                nant += 1
                index[i, 3] = nant  # Profundidade da barra na lateral
                if dados[i][2] == 0:  # Verifica se é barra terminal
                    lant += 1
                    nant = 0
                    nbf += 1
                    bf.append(dados[i][0])  # Determina vetor BF
                else:
                    lant = index[i, 1]  # Próxima barra está no mesmo nível
            else:  # Nível da barra já determinado
                lant = index[i, 1]
                nant = 1

            # Verifica se a barra deriva para outras laterais
            if len(dados[i]) > 3:
                for j in range(int(dados[i][3])):
                    j1 = 4 + j
                    derivada = dados[i][j1] - 1  # Ajuste de índice
                    index[derivada, 3] = 1  # Inicializa a profundidade da nova lateral ou sublateral
                    if j == 0:
                        index[derivada, 1] = index[i, 1] + j + 1  # Incrementa o nível da lateral
                        lnew = index[derivada, 1]
                        maux[lnew - 1] += 1  # Ajuste para índice Python
                        index[derivada, 2] = maux[lnew - 1]  # Profundidade da lateral

    # Determinação do vetor RBF
    rbf = bf[::-1]  # Inverte o vetor BF para obter RBF

    return index, rbf, nbf    


def run_fpd(topology_data, load_data, tol, max_it, ref,final_ref,v_inicial,ang_inicial):
    
    tol = np.float64(tol)
    max_it = int(max_it)
    ref = int(ref)
    final_ref = int(final_ref)
    ramos = set()

    # Leitura de dados de topologia, cargas e parâmetros da rede
    dados, carga, line, nbus, param = ler_dados(topology_data ,load_data, ref, final_ref)

    # Indexação das barras
    index, rbf, nbf = indexacao(dados, nbus)

    # Inicializações
    corr_barra = np.zeros(nbus,dtype=complex)
    corr_ramo  = np.zeros((nbus-1,nbus),dtype=complex)
    tensao = []

    # Inicializações adicionais
    pativ_ant = 0.0  # total de perdas ativas na iteração anterior
    preat_ant = 0.0  # total de perdas reativas na iteração anterior
    dpativ = 999.0   # diferença no valor da perda ativa entre duas iterações
    dpreat = 999.0   # diferença no valor da perda reativa entre duas iterações

    # Inicialização das tensões nas barras
    tensao = np.ones(nbus, dtype=complex) * (v_inicial * np.cos(ang_inicial * np.pi / 180) + 1j * v_inicial * np.sin(ang_inicial * np.pi / 180))
    iter = 0
    
    objetivo=[]
    # Processo iterativo
    while iter <= max_it and (dpativ >= tol or dpreat >= tol):
        iter += 1
        print(dpativ,dpreat)
        # Cálculo das correntes nas barras
        for i in range(nbus):
            corr_barra[i] = np.conj(complex(carga['ativ'][i], carga['reat'][i]) / tensao[i])

        # Somatórios acumulados (Backward)
        for i in range(nbf):
            bterm = rbf[i]
            bpara = bterm-1
            for j in range(index[bterm-1, 3]):
                if dados[bpara][1] != 0:
                    bde = dados[bpara][1]-1
                    corr_ramo[bde, bpara] = corr_barra[bpara]
                    corr_barra[bde] += corr_ramo[bde, bpara]
                    ramos.add((bde,bpara))
                    bpara = bde
                    
                    

        # Recalcula as tensões (Forward)
        perda_ativ = 0.0
        perda_reat = 0.0
        for i in range(nbus):
            bde = dados[i][1]-1  # barra anterior
            if bde > -1:
                bpara = i          # barra atual
                print(bde,bpara)
                if(bpara<bde): 
                    tensao[bpara] = tensao[bde] - corr_ramo[bpara, bde] * complex(line['r'][bpara, bde], line['x'][bpara, bde])
                    perda_ativ += line['r'][bpara, bde] * abs(corr_ramo[bpara, bde])**2
                    perda_reat += line['x'][bpara, bde] * abs(corr_ramo[bpara, bde])**2
                else:
                    tensao[bpara] = tensao[bde] - corr_ramo[bde, bpara] * complex(line['r'][bde, bpara], line['x'][bde, bpara])
                    perda_ativ += line['r'][bde, bpara] * abs(corr_ramo[bde, bpara])**2
                    perda_reat += line['x'][bde, bpara] * abs(corr_ramo[bde, bpara])**2

        # Determinação das diferenças entre as perdas calculadas em duas iterações subsequentes
        dpativ = abs(perda_ativ - pativ_ant)
        dpreat = abs(perda_reat - preat_ant)
        pativ_ant = perda_ativ
        preat_ant = perda_reat
        
        objetivo.append((dpativ,dpreat))
    
    tensao_dict = {'Barra':list(range(1,len(tensao)+1)),
                   'V':[np.abs(x) for x in tensao],
                   'Angulo (graus)': [np.angle(x,deg=True) for x in tensao]}
    
    tensao_pd = pd.DataFrame(tensao_dict)
    
    ramos = list(ramos)
    ramos.sort(key=lambda x:x[0])
    correntes = []
    for ramo in ramos :
        correntes.append(corr_ramo[ramo[0]][ramo[1]])
    corr_ramo_dict = {'Ramo' : [f'{x[0]+1} - {x[1]+1}' for x in ramos],
                     'Corrente (pu)' : [np.abs(x) for x in correntes],
                     'Ângulo (graus)' : [np.angle(x,deg=True) for x in correntes]}
    
    corr_ramo_pd = pd.DataFrame(corr_ramo_dict)
        
    return(tensao_pd, corr_ramo_pd, iter, param, carga, nbus,objetivo)