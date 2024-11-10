from multiple_injections import *
from scipy import integrate
from scipy import optimize
def score(T, Y, toxic=False):
    Pstar = Y[:,1:].sum(axis=1)
    if toxic:
        return integrate.trapezoid(Pstar, x=T)*Y[:,0].max() #required to take the min, otherwise the function will get the max C to 0 by injecting all at the end
    else:
        return integrate.trapezoid(Pstar, x=T)
    #return integrate.trapezoid(Y[:,1], x=T)

def times_to_injections(Temps_injections, C_i:float, deltas=False):
    if deltas:
        return [(t, C_i) for t in Temps_injections.cumsum()]
    else:
        return [(t, C_i) for t in Temps_injections]

def func_to_optimise(Temps_injections:np.array, args, deltas=True):
    C_tot, P, Q, Qp, t_stop, parameters, step, toxic, n_inj= args
    C_i = C_tot/n_inj
    injections = times_to_injections(Temps_injections, C_i, deltas)
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y, toxic)

def func_to_optimise_const(values:np.array, args):
    deltas=True
    Ti, delta = values
    C_tot, P, Q, Qp, t_stop, parameters, step, toxic, n_inj= args
    Temps_injections = np.array([Ti] + [delta]*(n_inj-1))
    C_i = C_tot/n_inj
    injections = times_to_injections(Temps_injections, C_i, deltas)
    T, Y = run_simulation(P, Q, Qp, injections, t_stop, parameters, step)
    return score(T, Y, toxic)

def find_optimal(n_inj, C_tot, P, Q, Qp, t_stop, Ns, toxic=False, const = False, parameters=parameters, step=0.2):
    """Returns Temps_injection as an array of the time deltas, not the absolute times"""
    args = C_tot, P, Q, Qp, t_stop, parameters, step, toxic, n_inj
    if const:
        bounds = (0, t_stop), (0, t_stop/n_inj)
        func = func_to_optimise_const
    else:
        bounds = tuple([(0, t_stop)] + [(0, t_stop/n_inj)]*(n_inj-1))
        func = func_to_optimise
    intiall_guess = optimize.brute(func, bounds, args=(args,), Ns=Ns, finish=None)
    try :
        result = optimize.minimize(func, intiall_guess, args=(args,), method='L-BFGS-B', bounds=bounds, options={"ftol":1e-15})
    except ValueError:
        raise ValueError(intiall_guess, bounds)
    if const:
        Ti, delta = result.x
        Temps_injections = np.array([Ti] + [delta]*(n_inj-1))
    else:
        Temps_injections = result.x
    s = func_to_optimise(Temps_injections, args, deltas=True)/t_stop
    return Temps_injections, s, result.message