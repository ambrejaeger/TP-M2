from model import *

def run_simulation(P0:float, Q0:float, Qp0:float, injections:list[tuple[float]], C0=0.0, parameters=parameters_PCV, step=0.2):
    t=0
    y = C, P, Q, Qp = C0, P0, Q0, Qp0
    Y = np.array([[C], [P], [Q], [Qp]])

    for delta_t, delta_C in injections:
        C += delta_C
        y = C, P, Q, Qp

        if delta_t > 0:
            Tstep = np.concatenate((np.arange(t, t+delta_t, step), t+delta_t))
            Ystep = np.concatenate(Y, odeint(derivees, Y[-1], T, args=(parameters_PCV,)))
            t += delta_t
        else :
            Tstep = np.array([t,])



    