import numpy as np

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

def Hohman_transfert():
    """Calcul des paramètres de transfert de Hohmann entre une orbite à 200km  et 300km."""
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

    "Partie 2: Deorbit maneuver "
    #Compute the propellant required to do a deorbit maneuver at 300km using a Hohmann transfer with a perigee at 10Km
    Ra = 300e3 + R_e
    Rp = 10e3 + R_e
    
    Dv_a2 = np.sqrt(K / Ra) * (1 - np.sqrt(2 /(1+Ra/Rp)))
    Dv = Dv_a
    m_final2 = M_tot / (np.exp( Dv / (I_sp * g) )) # Masse finale après consommation de carburant
    m_fuel2 = M_tot - m_final2  # Masse de carburant consommée

    #Compute the velocity and flight path angle of the capsule at re-entry assuming a re-entry point at 122km.
    h = 122e3
    re = h + R_e
    gamma_e = -6.2

    v_e = h / (re * np.cos(gamma_e))


    return Dv_a2, m_fuel2, m_final2, v_e

Dv_a2, m_fuel2, m_final2, v_e = Hohman_transfert()

print(Dv_a2, m_fuel2, m_final2, v_e)

