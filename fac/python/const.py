from math import *

Hartree_eV = 27.2113962            # Hartree in eV   
Rate_AU = 4.13413733E16            # Atomic Rate Unit in s-1
Rate_AU10 = 4.13413733E06          # Atomic Rate Unit in 10^10 s-1
Rate_AU12 = 4.13413733E04          # Atomic Rate Unit in 10^12 s-1
Area_AU20 = 2.80028560859E3        # Atomic Area Unit in 10^{-20} cm2
Alpha = 7.29735308E-3              # Fine Structure Constant
Ryd_eV = 13.6056981                # Rydberg in eV
RBohr = 0.529177249                # Bohr radius in A
FWHM = 2.35482005                  # conversion from sigma to FWHM
hc = 1.239842E4                    # hc in eV*A
hbc = 1.97327053E3                 # h_bar*c in eV*A
Me_eV = 5.1099906E5                # electron mass in eV
Me_keV = 5.1099906E2               # electron mass in keV
Mp_MeV = 9.38271998E2              # proton mass in MeV
Mp_keV = 9.38271998E5              # proton mass in keV
c  = 2.99792458E10                 # speed of light in cm/s
c10 = 2.99792458                   # speed of light in 10^10 cm/s
e = 1.60217733E-19                 # electron charge in Coulomb
e19 = 1.60217733                   # electron charge in 10^-19 Coulomb
erg_eV = 6.241506363E-13           # erg in eV
re = 2.81794092                    # electron classical radius in fm
sig_t = 6.6524616                  # Thompson cross section in 10^{-25} cm2
kb = 8.617385E-5                   # Boltzman constant in ev/K
Maxwellian = 1.12837967            # Constant before the maxwellian dist.

def EnergyH(z, n, l, j):
    if (j == 1):
        ka = (l + 1.0)
    elif (j == -1):
        ka = l
    else:
        raise 'j must be +1 or -1'

    if (n <= 0):
        raise 'n must be > 0'

    if (l < 0 or l >= n):
        raise 'l must be > 0 and < n'

    az = Alpha * z
    g = sqrt(ka*ka - az*az)
    np = n + g - ka
    x = az/np;
    e = (1.0/sqrt(1.0+x*x)) - 1.0
    e = e * Me_eV
    return e
