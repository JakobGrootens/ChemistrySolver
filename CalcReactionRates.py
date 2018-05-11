###########################
# Reaction rate functions #
###########################

#From Grackle Source (phys_constants.h)
#Boltzmann constant
kboltz = 1.38064852 * (10**-23)
#kboltz = 8.6173303 * 10**-5

#Mass of hydrogen
mh = 1.67262171 * (10**-24)
#Tiny
tiny =1.e-20

minimum_temperature = 1.0

import numpy as np
#  ---1:--       HI    + e   -> HII   + 2e
#  ---2:--       HII   + e   -> H     + p
#  ---3:--       HeI   + e   -> HeII  + 2e
#  ---4:--       HeII  + e   -> HeI   + p
#  ---5:--       HeII  + e   -> HeIII + 2e
#  ---6:--       HeIII + e   -> HeII  + p
#  ---7:--       HI    + e   -> HM    + p
#  ---8:--       HM    + HI  -> H2I*  + e
#  ---9:--       HI    + HII -> H2II  + p
#  ---10--       H2II  + HI  -> H2I*  + HII
#  ---11--       H2I   + HII -> H2II  + H
#  ---12--       H2I   + e   -> 2HI   + e
#  ---13--       H2I   + H   -> 3H
#  ---14--       HM    + e   -> HI    + 2e
#  ---15--       HM    + HI  -> 2H    + e
#  ---16--       HM    + HII -> 2HI
#  ---17--       HM    + HII -> H2II  + e
#  ---18--       H2II  + e   -> 2HI
#  ---19--       H2II  + HM  -> HI    + H2I
#  ---57--       HI    + HI  -> HII   + HI    + e
#  ---58--       HI    + HeI -> HII   + HeI   + e


def temperature(state, density, T):
    HI = state[0]; HM = state[1]; HII = state[2]; HeI = state[3]
    HeII = state[4]; HeIII = state[5]; H2I = state[6]; H2II = state[7]
    e = state[8]

    #Derived from Grackle source (calculate_temperature.c)
    temperature_units = 5000000 * (1.66 - 1) * mh / kboltz
    mu = (2*HI + 2*HM + 2*HII + 4*HeI + 4*HeII + 4*HeIII + 4*H2I + 4*H2II) / \
         ((2*HI + 2*HM + 2*HII + 4*HeI + 4*HeII + 4*HeIII + 4*H2I + 4*H2II) + e)
    return temperature_units * mu

#  ---1:--       HI    + e   -> HII   + 2e
def k1(T, T_eV, log_T_eV):
    rate = np.exp(-32.71396786375
          + 13.53655609057*log_T_eV
          - 5.739328757388*log_T_eV**2
          + 1.563154982022*log_T_eV**3
          - 0.2877056004391*log_T_eV**4
          + 0.03482559773736999*log_T_eV**5
          - 0.00263197617559*log_T_eV**6
          + 0.0001119543953861*log_T_eV**7
          - 2.039149852002e-6*log_T_eV**8)
    if(T_eV <= .8):
        rate = np.max([rate, tiny])
    return rate


#  ---2:--       HII   + e   -> H     + p
def k2(T):
    if(T < 1.0e9):
        rate = 4.881357e-6*T**(-1.5)* (1.+1.14813e2 * T**(-0.407))**(-2.242)
    else:
        rate = tiny
    return rate


#  ---3:--       HeI   + e   -> HeII  + 2e
def k3(T, T_eV, log_T_eV):
    if(T_eV > .8):
        rate = np.exp(-44.09864886561001
                 + 23.91596563469*log_T_eV
                 - 10.75323019821*log_T_eV**2
                 + 3.058038757198*log_T_eV**3
                 - 0.5685118909884001*log_T_eV**4
                 + 0.06795391233790001*log_T_eV**5
                 - 0.005009056101857001*log_T_eV**6
                 + 0.0002067236157507*log_T_eV**7
                 - 3.649161410833e-6*log_T_eV**8)
    else:
        rate = tiny
    return rate


#  ---4:--       HeII  + e   -> HeI   + p
def k4(T):
    rate = 1.26e-14 * (5.7067e5/T)**(0.75)
    return rate


#  ---5:--       HeII  + e   -> HeIII + 2e
def k5(T, T_eV, log_T_eV):
    if(T_eV > .8):
        rate = np.exp(-68.71040990212001
                 + 43.93347632635*log_T_eV
                 - 18.48066993568*log_T_eV**2
                 + 4.701626486759002*log_T_eV**3
                 - 0.7692466334492*log_T_eV**4
                 + 0.08113042097303*log_T_eV**5
                 - 0.005324020628287001*log_T_eV**6
                 + 0.0001975705312221*log_T_eV**7
                 - 3.165581065665e-6*log_T_eV**8)
    else:
        rate = tiny
    return rate


#  ---6:--       HeIII + e   -> HeII  + p
def k6(T):
    if(T < 1.0e9):
        rate = 7.8155e-5*T**(-1.5)* (1.+2.0189e2* T**(-0.407))**(-2.242)
    else:
        rate = tiny
    return rate


#  ---7:--       HI    + e   -> HM    + p
def k7(T):
    rate = 3.0e-16 * (T/3.e2)**0.95 * np.exp(-T/9.32e3)
    return rate


#  ---8:--       HM    + HI  -> H2I*  + e
def k8(T):
    rate = 1.35e-9*(T**9.8493e-2 + 3.2852e-1
            * T**5.5610e-1 + 2.771e-7 * T**2.1826) / (1. + 6.191e-3 * T**1.0461
            + 8.9712e-11 * T**3.0424
            + 3.2576e-14 * T**3.7741)
    return rate

#  ---bork bark:--       HI    + HII -> H2II  + p
def k9(T):
    if(T < 30.0):
        dog = 2.10e-20 * (T/30.0)**(-0.15)
    else:
        if(T > 3.2e4):
            tk9 = 3.2e4
        else:
            tk9 = T
        dog = 1e1**(-18.20 - 3.194 * np.log10(tk9) \
                + 1.786 * np.log10(tk9)**2 \
                - 0.2072 * np.log10(tk9) ** 3)
    return dog


#  ---10--       H2II  + HI  -> H2I*  + HII
def k10(T):
    rate = 6.0e-10
    return rate


#  ---11--       H2I   + HII -> H2II  + H
def k11(T, T_eV, log_T_eV):
    if(T_eV > .3):
        rate = np.exp(-24.24914687731536
                 + 3.400824447095291*log_T_eV
                 - 3.898003964650152*log_T_eV**2
                 + 2.045587822403071*log_T_eV**3
                 - 0.5416182856220388*log_T_eV**4
                 + 0.0841077503763412*log_T_eV**5
                 - 0.007879026154483455*log_T_eV**6
                 + 0.0004138398421504563*log_T_eV**7
                 - 9.36345888928611e-6*log_T_eV**8)
    else:
        rate = tiny
    return rate


#  ---12--       H2I   + e   -> 2HI   + e
def k12(T, T_eV):
    if(T_eV > .3):
        rate = 4.4886e-9*T**0.109127*np.exp(-101858./T)
    else:
        rate = tiny
    return rate


#  ---13--       H2I   + H   -> 3H
def k13(T, T_eV):
    if(T_eV > .3):
        rate = 1.0670825e-10*T_eV**2.012/(np.exp(4.463/T_eV)*(1.+0.2472* T_eV)**3.512)
    else:
        rate = tiny
    return rate


#  ---14--       HM    + e   -> HI    + 2e
def k14(T, T_eV, log_T_eV):
    if(T_eV > .04):
        rate = np.exp(-18.01849334273
                 + 2.360852208681*log_T_eV
                 - 0.2827443061704*log_T_eV**2
                 + 0.01623316639567*log_T_eV**3
                 - 0.03365012031362999*log_T_eV**4
                 + 0.01178329782711*log_T_eV**5
                 - 0.001656194699504*log_T_eV**6
                 + 0.0001068275202678*log_T_eV**7
                 - 2.631285809207e-6*log_T_eV**8)
    else:
        rate = tiny
    return rate


#  ---15--       HM    + HI  -> 2H    + e
def k15(T, T_eV, log_T_eV):
    if(T_eV > 0.1):
        rate = np.exp(-20.37260896533324
                     + 1.139449335841631*log_T_eV
                     - 0.1421013521554148*log_T_eV**2
                     + 0.00846445538663*log_T_eV**3
                     - 0.0014327641212992*log_T_eV**4
                     + 0.0002012250284791*log_T_eV**5
                     + 0.0000866396324309*log_T_eV**6
                     - 0.00002585009680264*log_T_eV**7
                     + 2.4555011970392e-6*log_T_eV**8
                     - 8.06838246118e-8*log_T_eV**9)
    else:
        rate = 2.56e-9*T_eV**1.78186
    return rate


#  ---16--       HM    + HII -> 2HI
def k16(T):
    rate = 2.4e-6*(1.+T/2e4)/np.sqrt(T)
    return rate


#  ---17--       HM    + HII -> H2II  + e
def k17(T):
     if (T > 1.0e4):
        rate = 4.0e-4*T**(-1.4)*np.exp(-15100./T)
     else:
        rate = 1.e-8*T**(-0.4)
     return rate


#  ---18--       H2II  + e   -> 2HI
def k18(T):
    if (T > 617):
        rate = 1.32e-6 * T ** (-0.76)
    else:
        rate = 1.e-8
    return rate


#  ---19--       H2II  + HM  -> HI    + H2I
def k19(T):
    rate = 5.e-7*np.sqrt(100./T)
    return rate


#  ---57--       HI    + HI  -> HII   + HI    + e
def k57(T):
    rate = 1.2e-17  * T**1.2e0 * np.exp(-1.578e5 / T)
    return rate

#  ---58--       HI    + HeI -> HII   + HeI   + e
def k58(T):
    rate = 1.75e-17 * T**1.3e0 * np.exp(-1.578e5 / T)
    return rate
