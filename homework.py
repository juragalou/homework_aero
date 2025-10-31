import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt




"""donné initial """

g = 9.81  # Accélération due à la gravité en m/s²
rho = 1.225  # Densité de l'air au niveau de la mer
R_e = 6371000  # Rayon de la Terre en mètres
C_d = 0.2  # Coefficient de traînée pour une sphère
M_module = 2460  # Masse du module lunaire en kg
M_tot = 4730  # Masse totale du vaisseau spatial en kg
T = 15.83e3  # Poussée du moteur en Newtons
I_sp = 266  # Impulsion spécifique en secondes
D = 2.3  # Diamètre du module en mètres
phi_i = 60 * (np.pi / 180)  # Angle initial en radians
R = 287.15 


def Hohman_transfert():

    "Partie 1: Space maneuver (service + reentry modules)"
    K = g * R_e**2  # Constante gravitationnelle terrestre

    Rp = R_e + 200e3  # Rayon de l'orbite initiale en mètres
    Ra = R_e + 300e3  # Rayon de l'orbite finale

    a = (Rp + Ra) / 2  # Demi-grand axe de l'orbite de transfert
    e = (Ra - Rp) / (Ra + Rp)  # Excentricité de l'orbite de transfert

    t = (np.pi * np.sqrt(a**3 / K))  # Temps de transfert en secondes

    Dv_p = np.sqrt(K / Rp) * (np.sqrt((2 * Ra / Rp)/(1+Ra/Rp)) - 1)  # Delta-v au périgée
    Dv_a = np.sqrt(K / Ra) * (1 - np.sqrt(2 /(1+Ra/Rp)))  # Delta-v à l'apogée

    Dv = Dv_p + Dv_a


    m_final = M_tot / (np.exp( Dv / (I_sp * g) )) # Masse finale après consommation de carburant
    m_fuel = M_tot - m_final  # Masse de carburant consommée
    m_debit = T / (I_sp * g)  # Débit massique du moteur

    omega_i = np.sqrt(K / Rp**3)  # Vitesse angulaire orbitale initiale
    omega_f = np.sqrt(K / Ra**3)  # Vitesse angulaire orbitale finale
    phi_f = np.pi - omega_i * t  # Angle final en radians

    T_wait = (phi_i - phi_f ) / (omega_i - omega_f)  # Temps d'attente avant l'injection

    t_impulsion = (m_fuel / m_debit)  # Durée de l'impulsion de poussée
    
    "on voit qu'on a un temps d'impulsion très court par rapport au temps de transfert, "

#-----------------------------------------

    "Partie 2: Deorbit maneuver "
    #Compute the propellant required to do a deorbit maneuver at 300km using a Hohmann transfer with a perigee at 10Km
    Ra2 = 300e3 + R_e
    Rp2 = 10e3 + R_e

    a2 = (Rp2 + Ra2) / 2 
    
    Dv_a2 = np.sqrt(K / Ra2) * (1 - np.sqrt(2 /(1+Ra2/Rp2))) #delta de desorbitation
    Dv = Dv_a2

    m_final2 = M_tot / (np.exp( Dv / (I_sp * g) )) # Masse finale après consommation de carburant
    m_fuel2 = M_tot - m_final2  # Masse de carburant consommée

    #Compute the velocity and flight path angle of the capsule at re-entry assuming a re-entry point at 122km.
    h = 122e3
    re = h + R_e
    #gamma_e = np.deg2rad(-6.2) #flight path angle  doit etre compris entre 5.5° et 6.9°



    va = np.sqrt(K * (2/Ra2 - 1/a2)) #vitesse à l'apogée 
    h_spec = Ra2 * va
    ve = np.sqrt(va**2 + 2*K*(1/re - 1/Ra2)) #vitesse au point d'entrée 
    gamma_e = -np.arccos(h_spec / (re * ve)) #flight path angle 
    gamma = np.rad2deg(gamma_e)






#---------------------------------------------------------

    "Partie 3: Flight path entry (reentry module only)"
    #Assuming a ballistic entry with the evolution of density in 3.a, and
    #using the initial conditions from question 2.b, plot the velocity
    #from reentry point with respect to altitude up to h=40km and the
    #deceleration (g loads) faced by the pilot. Mention in the report
    #values for 80, 60 and 40 km.

    h_40 = np.linspace(122e3, 40e3, 100)
    h_60 = np.linspace(122e3, 60e3, 100)
    h_80 = np.linspace(122e3, 80e3, 100)

    S = np.pi * (D/2)**2
    C_B = 2460 / (S * C_d) 

    rho_h = rho * np.exp(-0.1378e-3 * h_40)


    v_40= ve * np.exp((rho * (np.exp(-0.1378e-3* h_40)-np.exp(-0.1378e-3 * 122e3)))/(2*0.1378e-3*C_B*np.sin(gamma_e)))
    v_60 = ve * np.exp((rho * (np.exp(-0.1378e-3* h_60)-np.exp(-0.1378e-3 * 122e3)))/(2*0.1378e-3*C_B*np.sin(gamma_e)))
    v_80 = ve * np.exp((rho * (np.exp(-0.1378e-3* h_80)-np.exp(-0.1378e-3 * 122e3)))/(2*0.1378e-3*C_B*np.sin(gamma_e)))

    # F = m a --> F = trainée 

    Drag = 0.5 * rho_h * v_40 **2 * S * C_d

    acc = Drag / (2460 * g) 

    Q = 1.74e-4 * np.square((rho_h * v_40**3)/(D/2))

    I = simpson(-Q/(acc*g), v_40)#integre Q par rapport a V avec le changement de variable dt = dV/a(V)
    I1 = simpson(Q/(v_40 * np.sin(gamma_e)),h_40) #integre Q par rapport a h avec le changement de variable dt = dh/(V *sin(gamma_e))
    
    plt.figure(figsize=(8,6))
    plt.plot(h_80/1e3, v_80/1e3,   linewidth=8, label="jusqu'à 80 km")
    plt.plot(h_60/1e3, v_60/1e3, linewidth=4, label="jusqu'à 60 km")
    plt.plot(h_40/1e3, v_40/1e3, linewidth=2, label="jusqu'à 40 km")

    plt.gca().invert_xaxis()
    plt.xlabel("Altitude h (km)")
    plt.ylabel("Vitesse v (km/s)")
    plt.title("Rentrée balistique : vitesse en fonction de l’altitude")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8,6))
    plt.plot(h_40/1e3, acc, label= "acceleration")

    plt.gca().invert_xaxis()
    plt.xlabel("Altitude h (km)")
    plt.ylabel("g-loads")
    plt.title("Rentrée balistique : deceleration en fonction de l’altitude")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8,6))
    plt.plot(h_40/1e3, Q, label= "heat flux ")

    plt.gca().invert_xaxis()
    plt.xlabel("Altitude h (km)")
    plt.ylabel("heat flux")
    plt.title("Rentrée balistique : flux de chaleur en fonction de l’altitude")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()




    return I, I1




    



print(Hohman_transfert())

