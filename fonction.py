import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt



def density():

    g = 9.81       # m/s^2
    Rgas = 287.15  # J/(kg·K)

    # Bornes de couches (km -> m)
    hi = np.array([0.0, 11, 20, 32, 47, 51, 71, 84.85]) * 1000.0
    # Gradients (K/km -> K/m)
    L_hi = np.array([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]) * 1e-3
    # Températures de base (K) à h = hi[i]
    Ti = np.array([288.15, 216.650, 216.650, 228.650, 270.650, 270.650, 214.650, 186.946])
    # Pressions de base (Pa) à h = hi[i]
    p_base = np.array([101325, 22632, 5474.9, 868.02, 110.91, 66.939, 3.9564, 0.37338])

    def integrate_layer(i):
        h0, h1 = hi[i], hi[i+1]
        Li, Ti0 = L_hi[i], Ti[i]

        def dpdh(h, p):
            T = Ti0 + Li * (h - h0)
            return - (g * p) / (Rgas * T)

        h_eval = np.linspace(h0, h1, 300)
        sol = solve_ivp(dpdh, (h0, h1), [p_base[i]], t_eval=h_eval, rtol=1e-9, atol=1e-12)
        T = Ti0 + Li * (h_eval - h0)

        return h_eval, sol.y[0],T   # 1D arrays

    # Intégrer toutes les couches
    H_list, P_list, T_list = [], [], []
    for i in range(len(hi)-1):
        H_i, P_i,T_i = integrate_layer(i)
        H_list.append(H_i)
        P_list.append(P_i)
        T_list.append(T_i)


    H = np.concatenate(H_list)          # m
    P = np.concatenate(P_list)
    T = np.concatenate(T_list)         # Pa

    rho = P/ (Rgas * T )
    """
    # Tracé (h en km)
    plt.plot(H/1000.0, rho, label="rho(h) numérique")
    plt.xlabel("Altitude h (km)")
    plt.ylabel(r"Densité $\rho$ (kg/m$^3$)")
    plt.grid(True)
    plt.legend()
    plt.show()"""

    return H, rho
















def density_simple():

    g = 9.81       # m/s^2
    Rgas = 287.15  # J/(kg·K)
    rho_0 = 1.225

    h = np.linspace(0,84,100)
    rho = rho_0 * np.exp(-0.1378 * h)

   

    """
    # Tracé (h en km)
    plt.plot(h/1000.0, rho, label="rho(h) ")
    plt.xlabel("Altitude h (km)")
    plt.ylabel(r"Densité $\rho$ (kg/m$^3$)")
    plt.grid(True)
    plt.legend()
    plt.show()"""

    return h, rho


def compare_density():
    """Trace les deux courbes de densité sur un même graphique"""
    H1, rho1 = density()
    H2, rho2 = density_simple()

    plt.figure(figsize=(8, 5))
    plt.plot(H1/1000, rho1, label="Modèle ISA (réaliste)", linewidth=2)
    plt.plot(H2, rho2, "--", label="Modèle isotherme", linewidth=2)
    plt.xlabel("Altitude h (km)")
    plt.ylabel(r"Densité $\rho$ (kg/m$^3$)")
    plt.title("Comparaison des modèles de densité de l'air")
    plt.grid(True)
    plt.legend()
    plt.show()



print(compare_density())