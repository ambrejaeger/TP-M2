from multiple_injections import *
from scipy import integrate
from scipy import optimize
def score(T, Y):
    Pstar = Y[:,1:].sum(axis=1)
    return integrate.trapezoid(Pstar, x=T)
    #return integrate.trapezoid(Y[:,1], x=T)

def times_to_injections(Temps_injections:np.array, C_tot:float, deltas=False):
    if deltas:
        return [(t, 1/C_tot) for t in Temps_injections.cumsum()]
    else:
        return [(t, 1/C_tot) for t in Temps_injections]

def func_to_optimise(Temps_injections:np.array, args, deltas=True):
    C_tot, P, Q, Qp, t_stop, parameters, step= args
    injections = times_to_injections(Temps_injections, C_tot, deltas)
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y)

def find_optimal(n_inj, C_tot, P, Q, Qp, t_stop, maxfun, parameters=parameters, step=0.2):
    delta_min = 1.5
    bounds = tuple([(0, t_stop-delta_min*n_inj)] + [(delta_min, t_stop/n_inj)]*(n_inj-1))
    args = C_tot, P, Q, Qp, t_stop, parameters, step
    return optimize.direct(func_to_optimise, bounds, args=(args,), maxfun=maxfun)

if __name__ == "__main__":
    import time as t
    _, P0, Q0, Qp0 = y0_PCV
    n_inj = 2
    C_tot = 1
    t_stop = 120
    start = t.perf_counter()
    Result = find_optimal(n_inj, C_tot, P0, Q0, Qp0, t_stop, maxfun=10)
    
    stop = t.perf_counter()
    print("optimisation time", stop-start)
    Temps_injections = Result.x   
    print(Temps_injections)
    from matplotlib import cm
    X = np.linspace(0, 120, 20)
    Y = np.linspace(0, 120, 20)
    args = C_tot, P0, Q0, Qp0, t_stop, parameters, 0.2
    Z = np.array([[func_to_optimise(np.array([x, y]), args, deltas=False) for x in X] for y in Y])
    X, Y = np.meshgrid(X, Y)
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection='3d'))
    # Plot the surface.
    surf = ax2.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig2.colorbar(surf)

    plt.show()

