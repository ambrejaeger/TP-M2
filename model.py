
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Définir le système d'équations différentielles
y0_PCV = (1.0, 7.13, 41.2, 0.0) #Initial conditions
parameters = {
    "lambda_p":0.121,
    "k_P_Q":0.0295,
    "k_Qp_P":0.0031, 
    "delta_Qp":0.00867,
    "gamma":0.729, 
    "KDE":0.24, 
    "K":100
    }


def derivees(y:tuple, t, parameters:dict):
    """ Calcul des dérivées du système. 
    y : (C, P, Q, Qp)
    C : Concentration sanguine de l'agent
    P : population en prolifération
    Q : population en quiescence
    Qp : population en quiescence
    Pstar population totale

    """
    C, P, Q, Qp = y
    Pstar = P + Q + Qp
    dC = - parameters["KDE"] * C 
    dP = (parameters["lambda_p"] * P * (1-Pstar/parameters["K"]) #prolifération
        + parameters["k_Qp_P"]*Qp # Dormance endomagée -> prolifération
        - parameters["k_P_Q"]*P   # Induction de dormance
        - parameters["gamma"] * C* parameters["KDE"] * P # Mort de cellules
        )
    dQ = parameters["k_P_Q"]*P - parameters["gamma"]*C*parameters["KDE"]*Q
    dQp = parameters["gamma"]*C*parameters["KDE"]*Q - parameters["k_Qp_P"]*Qp - parameters["delta_Qp"]*Qp

    return dC, dP, dQ, dQp

def growth_no_treatment(y:tuple, t, parameters:tuple[float]):
    """ Calcul des dérivées du système. 
    y : (P, Q)
    C : Concentration sanguine de l'agent
    P : population en prolifération
    Q : population en quiescence
    Qp : population en quiescence
    Pstar population totale

    """
    C, P, Q, Qp = y
    lambda_p, k_P_Q, K = parameters

    Pstar = P + Q
    dC  = 0
    dP = - (lambda_p * P * (1-Pstar/K) #prolifération
        + k_P_Q*P   # Induction de dormance
        )
    dQ = - k_P_Q*P
    dQp = 0

    return dC, dP, dQ, dQp
 

# Time points to solve for
#t = np.linspace(0, 150, 100)

# Solve the system of equations
#y = odeint(derivees, y0_PCV, t, args=(parameters_PCV,))


def plot_MTD(derivees, y0, parameters, t):
    y = odeint(derivees, y0, t, args=(parameters,))
    plt.plot(t, y[:, 1] + y[:, 2] + y[:, 3], label='MTD')
    plt.xlabel('Time (months)')
    plt.ylabel('MTD (mm)')
    plt.legend()
    plt.show()
    return y





