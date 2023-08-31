import numpy as np


def p_Flow(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: tuple, n_bus: int) -> np.ndarray:
    gij = get_line_values('g', bus[0], bus[1], network_values)
    bij = get_line_values('b', bus[0], bus[1], network_values)
    gsi = get_line_values('gs', bus[0], bus[1], network_values)
    tetai = state_vector[bus[0] - 2] if bus[0] > 1 else 0
    tetaj = state_vector[bus[1] - 2] if bus[1] > 1 else 0
    tetaij = tetai - tetaj
    vi = state_vector[bus[0] + (n_bus - 2)]
    vj = state_vector[bus[1] + (n_bus - 2)]
    pFlow = vi*vi*(gsi + gij) - vi*vj*(gij*np.cos(tetaij) + bij*np.sin(tetaij))
    return pFlow

def q_Flow(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: tuple, n_bus: int) -> np.ndarray:
    gij = get_line_values('g', bus[0], bus[1], network_values)
    bij = get_line_values('b', bus[0], bus[1], network_values)
    bsi = get_line_values('bs', bus[0], bus[1], network_values)
    tetai = state_vector[bus[0] - 2] if bus[0] > 1 else 0
    tetaj = state_vector[bus[1] - 2] if bus[1] > 1 else 0
    tetaij = tetai - tetaj
    vi = state_vector[bus[0] + (n_bus - 2)]
    vj = state_vector[bus[1] + (n_bus - 2)]
    qFlow = -vi*vi*(bsi + bij) - vi*vj*(gij*np.sin(tetaij) - bij*np.cos(tetaij))
    return qFlow

def p_Inj(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: int, n_bus: int) -> np.ndarray:
    vi = state_vector[bus + (n_bus - 2)]
    tetai = state_vector[bus - 2] if bus > 1 else 0
    pInj = 0
    for j in range(n_bus):
        tetaj = state_vector[j - 1] if j > 0 else 0
        vj = state_vector[j + (n_bus - 1)]
        Gij = y_bar[bus-1, j].real
        Bij = y_bar[bus-1, j].imag
        tetaij = tetai - tetaj
        pInj += vi*vj*(Gij*np.cos(tetaij) + Bij*np.sin(tetaij))
    return pInj

def q_Inj(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: int, n_bus: int) -> np.ndarray:
    vi = state_vector[bus + (n_bus - 2)]
    tetai = state_vector[bus - 2] if bus > 1 else 0
    qInj = 0
    for j in range(n_bus):
        tetaj = state_vector[j - 1] if j > 0 else 0
        vj = state_vector[j + (n_bus - 1)]
        Gij = y_bar[bus-1, j].real
        Bij = y_bar[bus-1, j].imag
        tetaij = tetai - tetaj
        qInj += vi*vj*(Gij*np.sin(tetaij) - Bij*np.cos(tetaij))
    return qInj

def voltage(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: int, n_bus: int) -> np.ndarray:
    return state_vector[bus + (n_bus - 2)]



def get_line_values(parameter: str, bus_from: int, bus_to: int, network_values: dict) -> float:
    bus_from = int(bus_from - 1)
    bus_to = int(bus_to - 1)
    search_1 = (parameter, bus_from, bus_to)
    search_2 = (parameter, bus_to, bus_from)
    key = search_1 in network_values.keys()
    if key:
        return network_values.get(search_1)
    else:
        return network_values.get(search_2)


    
def derivada_p_Flow(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: tuple, n_bus: int) -> np.ndarray:
    state_amount = n_bus*2 -1
    pFlowVector = np.zeros(shape = (state_amount, ))
    for i in range(state_amount):
        if i < (n_bus-1):
            if i + 2 in bus:
                vi = state_vector[bus[0] + n_bus -2]
                vj = state_vector[bus[1] + n_bus -2]
                gij = get_line_values('g', bus[0], bus[1], network_values)
                bij = get_line_values('b', bus[0], bus[1], network_values)
                tetai = state_vector[bus[0] - 2] if bus[0] > 1 else 0
                tetaj = state_vector[bus[1] - 2] if bus[1] > 1 else 0
                pFlowVector[i] = vi*vj*(gij*np.sin(tetai-tetaj) - bij*np.cos(tetai-tetaj)) if i + 2 == bus[0] else -vi*vj*(gij*np.sin(tetai-tetaj) - bij*np.cos(tetai-tetaj))
            else:
                pFlowVector[i] = 0
        else:
            actual_bus = i - (n_bus - 2)
            if actual_bus in bus:
                vi = state_vector[bus[0] + n_bus -2]
                vj = state_vector[bus[1] + n_bus -2]
                gij = get_line_values('g', bus[0], bus[1], network_values)
                bij = get_line_values('b', bus[0], bus[1], network_values)
                gsi = get_line_values('gs', bus[0], bus[1], network_values)
                tetai = state_vector[bus[0] - 2] if bus[0] > 1 else 0
                tetaj = state_vector[bus[1] - 2] if bus[1] > 1 else 0
                pFlowVector[i] = -vj*(gij*np.cos(tetai-tetaj) + bij*np.sin(tetai-tetaj)) + 2*(gij + gsi)*vi if actual_bus == bus[0] else -vi*(gij*np.cos(tetai-tetaj) + bij*np.sin(tetai-tetaj))
            else:
                pFlowVector[i] = 0
    return pFlowVector

def derivada_q_Flow(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: tuple, n_bus: int) -> np.ndarray:
    state_amount = n_bus*2 -1
    qFlowVector = np.zeros(shape = (state_amount, ))
    for i in range(state_amount):
        if i < (n_bus-1):
            if i + 2 in bus:
                vi = state_vector[bus[0] + n_bus -2]
                vj = state_vector[bus[1] + n_bus -2]
                gij = get_line_values('g', bus[0], bus[1], network_values)
                bij = get_line_values('b', bus[0], bus[1], network_values)
                tetai = state_vector[bus[0] - 2] if bus[0] > 1 else 0
                tetaj = state_vector[bus[1] - 2] if bus[1] > 1 else 0
                qFlowVector[i] = -vi*vj*(gij*np.cos(tetai-tetaj) + bij*np.sin(tetai-tetaj)) if i + 2 == bus[0] else vi*vj*(gij*np.cos(tetai-tetaj) + bij*np.sin(tetai-tetaj))
            else:
                qFlowVector[i] = 0
        else:
            actual_bus = i - (n_bus - 2)
            if actual_bus in bus:
                vi = state_vector[bus[0] + n_bus -2]
                vj = state_vector[bus[1] + n_bus -2]
                gij = get_line_values('g', bus[0], bus[1], network_values)
                bij = get_line_values('b', bus[0], bus[1], network_values)
                bsi = get_line_values('bs', bus[0], bus[1], network_values)
                tetai = state_vector[bus[0] - 2] if bus[0] > 1 else 0
                tetaj = state_vector[bus[1] - 2] if bus[1] > 1 else 0
                qFlowVector[i] = -vj*(gij*np.sin(tetai-tetaj) - bij*np.cos(tetai-tetaj)) - 2*(bij + bsi)*vi if actual_bus == bus[0] else -vi*(gij*np.sin(tetai-tetaj) - bij*np.cos(tetai-tetaj))
            else:
                qFlowVector[i] = 0
    return qFlowVector

def derivada_p_Inj(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: int, n_bus: int) -> np.ndarray:
    state_amount = n_bus*2 -1
    pInjVector = np.zeros(shape = (state_amount, ))
    for i in range(state_amount):
        if i < (n_bus-1):
            if i + 2 == bus:
                vi = state_vector[bus + n_bus -2]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Bii = y_bar[bus -1, bus - 1].imag
                bus_looking = range(n_bus)
                aux_p = 0
                for j in bus_looking:
                    vj = state_vector[j + n_bus - 1]
                    Gij = y_bar[i+1, j].real
                    Bij = y_bar[i+1, j].imag
                    tetaj = state_vector[j - 1] if j > 0 else 0
                    tetaij = tetai - tetaj
                    aux_p += vi*vj*(Bij*np.cos(tetaij) - Gij*np.sin(tetaij))
                pInjVector[i] = aux_p - vi*vi*Bii
            else:
                vi = state_vector[bus + n_bus -2]
                vj = state_vector[i + n_bus]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Gij = y_bar[bus - 1, i+1].real
                Bij = y_bar[bus - 1, i+1].imag
                tetaj = state_vector[i]
                tetaij = tetai - tetaj
                pInjVector[i] = vi*vj*(Gij*np.sin(tetaij) - Bij*np.cos(tetaij))
        else:
            actual_bus = i - (n_bus - 2)
            if actual_bus == bus:
                vi = state_vector[bus + n_bus -2]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Gii = y_bar[bus -1, bus - 1].real
                bus_looking = range(n_bus)
                aux_p = 0
                for j in bus_looking:
                    vj = state_vector[j + n_bus - 1]
                    Gij = y_bar[actual_bus-1, j].real
                    Bij = y_bar[actual_bus-1, j].imag
                    tetaj = state_vector[j - 1] if j > 0 else 0
                    tetaij = tetai - tetaj
                    aux_p += vj*(Bij*np.sin(tetaij) + Gij*np.cos(tetaij))
                pInjVector[i] = aux_p + vi*Gii
            else:
                vi = state_vector[bus + n_bus -2]
                vj = state_vector[i]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Gij = y_bar[bus - 1, actual_bus-1].real
                Bij = y_bar[bus - 1, actual_bus-1].imag
                tetaj = state_vector[actual_bus-2] if actual_bus > 1 else 0
                tetaij = tetai - tetaj
                pInjVector[i] = vi*(Gij*np.cos(tetaij) + Bij*np.sin(tetaij))
    return pInjVector

def derivada_q_Inj(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: int, n_bus: int) -> np.ndarray:
    state_amount = n_bus*2 -1
    qInjVector = np.zeros(shape = (state_amount, ))
    for i in range(state_amount):
        if i < (n_bus-1):
            if i + 2 == bus:
                vi = state_vector[bus + n_bus -2]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Gii = y_bar[bus -1, bus - 1].real
                bus_looking = range(n_bus)
                aux_p = 0
                for j in bus_looking:
                    vj = state_vector[j + n_bus - 1]
                    Gij = y_bar[i+1, j].real
                    Bij = y_bar[i+1, j].imag
                    tetaj = state_vector[j - 1] if j > 0 else 0
                    tetaij = tetai - tetaj
                    aux_p += vi*vj*(Bij*np.sin(tetaij) + Gij*np.cos(tetaij))
                qInjVector[i] = aux_p - vi*vi*Gii
            else:
                vi = state_vector[bus + n_bus -2]
                vj = state_vector[i + n_bus]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Gij = y_bar[bus - 1, i+1].real
                Bij = y_bar[bus - 1, i+1].imag
                tetaj = state_vector[i]
                tetaij = tetai - tetaj
                qInjVector[i] = -vi*vj*(Gij*np.cos(tetaij) + Bij*np.sin(tetaij))
        else:
            actual_bus = i - (n_bus - 2)
            if actual_bus == bus:
                vi = state_vector[bus + n_bus -2]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Bii = y_bar[bus -1, bus - 1].imag
                bus_looking = range(n_bus)
                aux_p = 0
                for j in bus_looking:
                    vj = state_vector[j + n_bus - 1]
                    Gij = y_bar[actual_bus-1, j].real
                    Bij = y_bar[actual_bus-1, j].imag
                    tetaj = state_vector[j - 1] if j > 0 else 0
                    tetaij = tetai - tetaj
                    aux_p += vj*(Gij*np.sin(tetaij) - Bij*np.cos(tetaij))
                qInjVector[i] = aux_p - vi*Bii
            else:
                vi = state_vector[bus + n_bus -2]
                vj = state_vector[i]
                tetai = state_vector[bus - 2] if bus > 1 else 0
                Gij = y_bar[bus - 1, actual_bus-1].real
                Bij = y_bar[bus - 1, actual_bus-1].imag
                tetaj = state_vector[actual_bus-2] if actual_bus > 1 else 0
                tetaij = tetai - tetaj
                qInjVector[i] = vi*(Gij*np.sin(tetaij) - Bij*np.cos(tetaij))
    return qInjVector

def derivada_voltage(y_bar: np.ndarray, network_values: dict, state_vector: np.ndarray, bus: int, n_bus: int) -> np.ndarray:
    state_amount = n_bus*2 -1
    vVector = np.zeros(shape = (state_amount, ))
    for i in range(state_amount):
        if i < (n_bus-1):
            vVector[i] = 0
        else:
            actual_bus = i - (n_bus - 2)
            if actual_bus == bus:
                vVector[i] = 1
            else:
                vVector[i] = 0
    return vVector

