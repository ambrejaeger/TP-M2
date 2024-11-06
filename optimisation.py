from multiple_injections import *
from scipy import integrate
from scipy.optimize import brute
def score(T, Y):
    Pstar = Y[:,1:].sum(axis=1)
    return integrate.trapezoid(Pstar, x=T)

def func_to_optimise(Temps_injections:np.array, args):
    C_tot, P, Q, Qp, t_stop, parameters, step= args
    injections = [(t, 1/C_tot) for t in Temps_injections]
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y)

def find_optimal(n_inj, C_tot, fin_tmt, P, Q, Qp, t_stop, parameters=parameters_PCV, step=0.2):
    ranges = tuple([(0, fin_tmt)]*n_inj)
    args = C_tot, P, Q, Qp, t_stop, parameters, step
    return brute(func_to_optimise, ranges, args=(args,), Ns=20, full_output=True, finish=None)

if __name__ == "__main__":
    import time as t
    injections = [(10, 1), (13, 1), (16, 1), (19, 1), (21,1)]
    _, P0, Q0, Qp0 = y0_PCV
    start = t.perf_counter()
    Temps_injections, s, grid, Jout = find_optimal(3, 1, 120, P0, Q0, Qp0, 120)
    stop = t.perf_counter()
    print("optimisation time", stop-start)
    print(Temps_injections, s)
    plt.show()