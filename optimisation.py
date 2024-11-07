from multiple_injections import *
from scipy import integrate
from scipy import optimize
def score(T, Y, toxic=False):
    Pstar = Y[:,1:].sum(axis=1)/T[-1]
    if toxic:
        return integrate.trapezoid(Pstar, x=T)/Y[:,0].max()
    else:
        return integrate.trapezoid(Pstar, x=T)
    #return integrate.trapezoid(Y[:,1], x=T)

def times_to_injections(Temps_injections:np.array, C_i:float, deltas=False):
    if deltas:
        return [(t, C_i) for t in Temps_injections.cumsum()]
    else:
        return [(t, C_i) for t in Temps_injections]

def func_to_optimise(Temps_injections:np.array, args, deltas=True):
    C_tot, P, Q, Qp, t_stop, parameters, step, toxic= args
    C_i = C_tot/len(Temps_injections)
    injections = times_to_injections(Temps_injections, C_i, deltas)
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y, toxic)

def find_optimal(n_inj, C_tot, P, Q, Qp, t_stop, Ns, toxic=False, parameters=parameters, step=0.2):
    delta_min = 0.0
    bounds = tuple([(0, t_stop-delta_min*n_inj)] + [(delta_min, t_stop/n_inj)]*(n_inj-1))
    args = C_tot, P, Q, Qp, t_stop, parameters, step, toxic
    return optimize.brute(func_to_optimise, bounds, args=(args,), Ns=Ns, full_output=True)

def func_to_optimise_const(Temps_injections:np.array, args):
    deltas=True
    Ti, delta = Temps_injections
    C_tot, P, Q, Qp, t_stop, parameters, step, n_inj= args
    C_i = C_tot/n_inj
    injections = times_to_injections(Temps_injections, C_i, deltas)
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y)

def find_optimal_const(n_inj, C_tot, P, Q, Qp, t_stop, toxic, Ns, parameters=parameters, step=0.2):
    bounds = tuple((0, t_stop), (0, t_stop/n_inj))
    args = C_tot, P, Q, Qp, t_stop, parameters, step, toxic, n_inj
    return optimize.brute(func_to_optimise, bounds, args=(args,), Ns=Ns, full_output=True)

if __name__ == "__main__":
    import time as t
    _, P0, Q0, Qp0 = y0_PCV
    n_inj = 3
    C_tot = 1
    t_stop = 120
    step = 0.2
    args = C_tot, P0, Q0, Qp0, t_stop, parameters, step
    for Ns in range(0,0):
        print(Ns)
        for n_inj in range(1, 5):
            C_tot = n_inj
            print("n :", n_inj)
            start = t.perf_counter()
            Temps_injections, s, _, _ = find_optimal(n_inj, C_tot, P0, Q0, Qp0, t_stop, Ns=Ns, step=0.02)
            stop = t.perf_counter()
            print("Score : ", s)
            Temps_injections = Temps_injections.cumsum()
            #s = func_to_optimise(Temps_injections, args)
            print(Temps_injections)
            print("optimisation time :", stop-start)




plot_score_space = True
if plot_score_space:
    from matplotlib import cm
    X = np.linspace(0, 60, 20)
    Y = np.linspace(0, 60, 20)
    args = 1, P0, Q0, Qp0, t_stop, parameters, 0.2
    Z = np.array([[func_to_optimise(np.array([x, y]), args, deltas=False) for x in X] for y in Y])
    X, Y = np.meshgrid(X, Y)
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection='3d'))
    # Plot the surface
    surf = ax2.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
    fig2.colorbar(surf)

    plt.show()