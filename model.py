#Test comment Paul
# Définir le système d'équations différentielles
y0_PCV = (1.0, 7.13, 41.2, 0.0)
parameters_PCV = (0.121, 0.0295, 0.0031, 0.00867, 0.729, 0.729, 0.24, 100)

def derivees(y:tuple, parameters:tuple[float]):
    """ Calcul des dérivées partielles. 
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

if __name__=="__main__":

    dy = derivees(y0_PCV, parameters_PCV)
    print(dy)