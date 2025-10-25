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

    omega_i = np.sqrt(K / Rp**3)  # Vitesse angulaire orbitale initiale
    omega_f = np.sqrt(K / Ra**3)  # Vitesse angulaire orbitale finale
    phi_f = np.pi - omega_i * t  # Angle final en radians

    T_wait = (phi_i - phi_f ) / (omega_i - omega_f)  # Temps d'attente avant l'injection


    return t,a ,e, Dv_p, Dv_a, Dv_p + Dv_a, m_fuel, T_wait

t, a, e, Dv_p, Dv_a, Dv_tot, m_fuel, T_wait = Hohman_transfert()

print(f"Temps de transfert: {t/60:.2f} minutes")
print(f"Demi-grand axe: {a/1e3:.2f} km")
print(f"Excentricité: {e:.4f}")
print(f"Delta-v au périgée: {Dv_p:.2f} m/s")
print(f"Delta-v à l'apogée: {Dv_a:.2f} m/s")
print(f"Delta-v total: {Dv_tot:.2f} m/s")
print(f"Masse de carburant consommée: {m_fuel:.2f} kg")
print(f"Temps d'attente avant l'injection: {T_wait/60**2:.2f} heures")