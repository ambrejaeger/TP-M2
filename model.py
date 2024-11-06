
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Définir le système d'équations différentielles


def derivees(y:tuple, t, parameters:tuple[float]):
    """ Calcul des dérivées du système. 
    y : (C, P, Q, Qp)
    C : Concentration sanguine de l'agent
    P : population en prolifération
    Q : population en quiescence
    Qp : population en quiescence
    Pstar population totale

    """
    C, P, Q, Qp = y
    lambda_p, k_P_Q, k_Qp_P, delta_Qp, gamma_P, gamma_Q, KDE, K = parameters

    Pstar = P + Q + Qp
    dC = - KDE * C 
    dP = (lambda_p * P * (1-Pstar/K) #prolifération
        + k_Qp_P*Qp # Dormance endomagée -> prolifération
        - k_P_Q*P   # Induction de dormance
        - gamma_P * C* KDE * P # Mort de cellules
        )
    dQ = k_P_Q*P - gamma_Q*C*KDE*Q
    dQp = gamma_Q*C*KDE*Q - k_Qp_P*Qp - delta_Qp*Qp

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



plot_MTD(derivees, y0_PCV, parameters_PCV, np.linspace(0, 100, 1000))


