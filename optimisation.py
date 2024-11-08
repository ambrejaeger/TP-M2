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

if __name__ == "__main__":
    import time as t
    _, P0, Q0, Qp0 = y0_PCV
    n_inj = 6
    C_tot = 1
    t_stop = 40
    step = 0.2
    toxic = True
    Proto_s = []
    Single_s = []
    Const_s = []
    Optimal_s = []

    Const_ti = []
    Const_delt = []
    n_max = 7
    opti_times = np.empty((n_max,n_max,))
    opti_times.fill(np.nan)

    const_times = np.empty((n_max,n_max,))
    const_times.fill(np.nan)
    for n_inj in range(1, n_max+1):
        C_tot = n_inj
        args = C_tot, P0, Q0, Qp0, t_stop, parameters, step, toxic, n_inj
        
        print("n :", n_inj)
        start = t.perf_counter()

        
        
        #source protocol
        Temps_injections_proto = np.array([0] + [1.5]*(n_inj-1))
        s_proto = func_to_optimise(Temps_injections_proto, args, deltas=True)/t_stop
        print("Protocol :".ljust(12), Temps_injections_proto, s_proto)
        Proto_s.append(s_proto)

        #Single optimal injection
        Temps_injections_unique, s_unique, message = find_optimal(1, C_tot, P0, Q0, Qp0, t_stop, Ns=300, toxic=toxic, const=False, step=0.02)
        print("Single :".ljust(12), Temps_injections_unique, s_unique, message)
        Single_s.append(s_unique)
        if n_inj>1:
            #constant delta
            Temps_injections_const, s_const, message = find_optimal(n_inj, C_tot, P0, Q0, Qp0, t_stop, Ns=50, toxic=toxic, const=True,  step=0.02)
            print("Constant :".ljust(12), Temps_injections_const, s_const, message)
            
            Const_ti.append(Temps_injections_const[0])
            Const_delt.append(Temps_injections_const[1])
            
            if n_inj>2:
                #absolute_optimum
                bounds = tuple([(0, t_stop)] + [(0, t_stop/n_inj)]*(n_inj-1))
                Result = optimize.minimize(func_to_optimise, Temps_injections_const, args=(args,), method='L-BFGS-B', bounds=bounds, options={"ftol":1e-15})
                Temps_injections_opt = Result.x
                s_opt = func_to_optimise(Temps_injections_opt, args, deltas=True)/t_stop
                print("Optimal :".ljust(12), Temps_injections_opt, s_opt, Result.message)                
            else:
                Temps_injections_opt = Temps_injections_const
                s_opt = s_const
        else:
            Temps_injections_opt = Temps_injections_const = Temps_injections_unique
            s_opt = s_const = s_unique
        const_times[n_inj-1,:n_inj] = Temps_injections_const.cumsum()
        opti_times[n_inj-1,:n_inj] = Temps_injections_opt.cumsum()
        Const_s.append(s_const)
        Optimal_s.append(s_opt)



        stop = t.perf_counter()
        
        print("optimisation time :", stop-start)

    fig1, ax1 = plt.subplots()
    X = list(range(1, n_max+1))
    fig1.suptitle("Temps total de simulation "+ str(t_stop))
    for data, label in [(Proto_s, "Protocole clinique"), (Single_s, "Une injection equivalente"), (Const_s, "Frequence constante"), (Optimal_s, "optimal")]:
        print(label, len(data), len(X))
        ax1.plot(X, data, label=label)
    ax1.set_xlabel("Nombre d'injections totales")
    ax1.set_ylabel("Mean_diameter*Cmax")
    ax1.legend()

    fig2, ax2 = plt.subplots()
    ax2.plot(X[1:], Const_ti, label = "First injection")
    ax2.plot(X[1:], Const_delt, label = "Period")
    ax2.set_xlabel("n_inj")
    ax2.set_ylabel("time")
    ax2.legend()

    fig3, ax3 = plt.subplots()
    ax3.plot(X, opti_times, color="b", label="optimal")
    ax3.plot(X, const_times, color="r", label="constant")
    ax3.set_xlabel("n_inj")
    ax3.set_ylabel("time")




    plot_score_space = False
    if plot_score_space:
        from matplotlib import cm
        X = np.linspace(0, 60, 20)
        Y = np.linspace(0, 60, 20)
        args = 2, P0, Q0, Qp0, t_stop, parameters, 0.2, False, 2
        Z = np.array([[func_to_optimise(np.array([x, y]), args, deltas=False) for x in X] for y in Y])/120
        X, Y = np.meshgrid(X, Y)
        fig2, ax2 = plt.subplots(subplot_kw=dict(projection='3d'))
        # Plot the surface
        surf = ax2.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)
        fig2.colorbar(surf)

    plt.show()