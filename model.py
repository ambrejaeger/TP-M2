#Test comment Paul
# Définir le système d'équations différentielles



def derivees(y:tuple, parameters:tuple[float]):
    """ Calcul des dérivées partielles. 
    y : (C, P, Q, Qp)
    C : Concentration sanguine de l'agent
    P : population en prolifération
    Q : popualtion en quiescence
    Qp : population en quiescence
    Pstar population totale

    """
    C, P, Q, Qp = y
    KDE, lambda_p, K, k_Qp_P, k_P_Q, gamma_P, gamma_Q = parameters

    Pstar = P + Q + Qp
    dC = - KDE * C 
    dP = (lambda_p * P * (1-Pstar/K) #prolifération
        + k_Qp_P*Qp # Dormance endomagée -> prolifération
        - k_P_Q*P   # Induction de dormance
        - gamma_P * C* KDE * P # Mort de cellules
        )
    dQ = k_P_Q*P - gamma_Q*C*KDE*Q
    dQp = gamma_Q*C*KDE*Q

    return dC, dP, dQ, dQp