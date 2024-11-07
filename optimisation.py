from multiple_injections import *
from scipy import integrate
from scipy import optimize
def score(T, Y):
    Pstar = Y[:,1:].sum(axis=1)
    return integrate.trapezoid(Pstar, x=T)
    #return integrate.trapezoid(Y[:,1], x=T)

def times_to_injections(Temps_injections:np.array, C_i:float, deltas=False):
    if deltas:
        return [(t, C_i) for t in Temps_injections.cumsum()]
    else:
        return [(t, C_i) for t in Temps_injections]

def func_to_optimise(Temps_injections:np.array, args, deltas=True):
    C_tot, P, Q, Qp, t_stop, parameters, step= args
    C_i = C_tot/len(Temps_injections)
    injections = times_to_injections(Temps_injections, C_i, deltas)
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y)

def find_optimal(n_inj, C_tot, P, Q, Qp, t_stop, Ns, parameters=parameters, step=0.2):
    delta_min = 0.0
    bounds = tuple([(0, t_stop-delta_min*n_inj)] + [(delta_min, t_stop/n_inj)]*(n_inj-1))
    args = C_tot, P, Q, Qp, t_stop, parameters, step
    return optimize.brute(func_to_optimise, bounds, args=(args,), Ns=Ns, full_output=True)

if __name__ == "__main__":
    import time as t
    _, P0, Q0, Qp0 = y0_PCV
    n_inj = 3
    C_tot = 1
    t_stop = 70
    step = 0.2
    args = C_tot, P0, Q0, Qp0, t_stop, parameters, step
    for n_inj in range(1, 2):
        print("n :", n_inj)
        start = t.perf_counter()
        Temps_injections, s, grid, SGrid = find_optimal(n_inj, C_tot, P0, Q0, Qp0, t_stop, Ns=2000, step=0.02)
        stop = t.perf_counter()
        print("Score : ", s)
        #Temps_injections = Temps_injections.cumsum()
        injections = times_to_injections(Temps_injections, C_tot/n_inj, deltas=True)
        T, Y = run_simulation(P0, Q0, Qp0, injections, 120)
        plot_simulation(T, Y)
        s = score(T, Y)
        plt.title(str(s))
        #s = func_to_optimise(Temps_injections, args)
        print(Temps_injections)
        print("optimisation time :", stop-start)
        print("Score : ", s)




plot_score_space = False
if plot_score_space:
    from matplotlib import cm
    X = np.linspace(15, 35, 20)
    Y = np.linspace(15, 35, 20)
    args = C_tot, P0, Q0, Qp0, t_stop, parameters, 0.2
    Z = np.array([[func_to_optimise(np.array([x, y]), args, deltas=False) for x in X] for y in Y])
    X, Y = np.meshgrid(X, Y)
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection='3d'))
    # Plot the surface.
    surf = ax2.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig2.colorbar(surf)

    plt.show()