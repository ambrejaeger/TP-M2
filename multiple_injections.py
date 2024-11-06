from model import *

def run_simulation(P:float, Q:float, Qp:float, injections:list[tuple[float]], t_stop, parameters=parameters_PCV, step=0.2):
    """
    Takes initial parameters of the model of the tumor (P, Q and Qp), then list of injections
    List of injections to be formated as a tuple of floats: 
    the first representing the absolute time of injection,
    the second the concentration injected.

    Initial concentration is 0.

    Injection that happen after t_stop will be disregarded as if they never happened
    """
    events = injections.copy() + [(t_stop, None)] 
    #copying the list and adding the stop condition
    #copying is required to prevent the list being modified in whatever function created it
    #None will be recognised as a signal to exit the loop
    events.sort(key=lambda event: event[0])

    old_t=0
    C = 0
    #initialisation des valeurs de sorties
    Y = np.array([[C, P, Q, Qp]])
    T = np.array([0.0])

    for new_t, delta_C in events:
        if new_t > old_t: #there is a time-span that needs to be simulated
            y = C, P, Q, Qp
            n_points = int((new_t - old_t)/step) + 2
            Tstep = np.linspace(old_t, new_t, n_points, np.array([new_t]))
            #arange ne comprends pas la valeur de temps finale, on doit donc l'ajouter manuellement afin qu'il soit compris
            Ystep = odeint(derivees, y, Tstep, args=(parameters,))
            C, P, Q, Qp = Ystep[-1] #mise à jour des populations
            
            #on ajoute les valeurs qui viennent d'être calculées à nos variables d'export
            T = np.concatenate((T, Tstep))
            Y = np.concatenate((Y, Ystep), axis=0)
            old_t = new_t

        if delta_C is None: #we have reached t_stop
            return T, Y
        

        C += delta_C      

def plot_simulation(t, y): #je l'ai remise, au moins pour voir ce que donne une simulation avec injections multiples
    Pstar = y[:, 1] + y[:, 2] + y[:, 3]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    ax1.plot(t, Pstar, label='Total tumor size')
    ax1.plot(t, y[:, 1], ls="--", label='Proliferative')
    ax1.plot(t, y[:, 2], ls="--", label='Dormant')
    ax1.plot(t, y[:, 3], ls="--", label='Damaged dormant')
    ax1.set_xlabel("time (mo)")
    ax1.set_ylabel("Tumor size (mm)")

    ax2.plot(t, y[:, 0], "k")
    ax2.set_ylabel("Drug concentration (AU)")
    fig.legend()
    return fig

if __name__ == "__main__":
    injections = [(10, 1), (13, 1), (16, 1), (19, 1), (21,1)]
    _, P0, Q0, Qp0 = y0_PCV
    T, Y = run_simulation(P0, Q0, Qp0, injections, 100)
    fig = plot_simulation(T, Y)
    plt.show()