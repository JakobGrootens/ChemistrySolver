import CalcReactionRates as rr
import scipy.integrate as sint
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import time
import sys

T = 0
energy = 0

#Calculates and returns the reaction rates from a given Temperature
def calculate_reactionrates(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)

    k1 = rr.k1(T, T_eV, log_T_eV)
    k2 = rr.k2(T)
    k3 = rr.k3(T, T_eV, log_T_eV)
    k4 = rr.k4(T)
    k5 = rr.k5(T, T_eV, log_T_eV)
    k6 = rr.k6(T)
    k7 = rr.k7(T)
    k8 = rr.k8(T)
    k9 = rr.k9(T)
    k10 = rr.k10(T)
    k11 = rr.k11(T, T_eV, log_T_eV)
    k12 = rr.k12(T, T_eV)
    k13 = rr.k13(T, T_eV)
    k14 = rr.k14(T, T_eV, log_T_eV)
    k15 = rr.k15(T, T_eV, log_T_eV)
    k16 = rr.k16(T)
    k17 = rr.k17(T)
    k18 = rr.k18(T)
    k19 = rr.k19(T)
    k21 = rr.k21(T)
    k22 = rr.k22(T)
    k57 = rr.k57(T)
    k58 = rr.k58(T)

    return k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, \
    k11, k12, k13, k14, k15, k16, k17, k18,  \
    k21, k22, k19, k57, k58


#Calculates the change in species for each slice
def updater(t, state):

    HI = state[0]; HM = state[1]; HII = state[2]; HeI = state[3]
    HeII = state[4]; HeIII = state[5]; H2I = state[6]; H2II = state[7]
    e = HII + HeII + 2*HeIII + H2II - HM

    #Expressions made from RHSGenerator.py
    dHIdt =  -2*H2I*HI*k21 + H2I*HII*k11 + H2I*e*k12 - H2II*HI*k10 + H2II*HM*k19 + H2II*e*k18 - \
             HI*HII*k9 - HI*HM*k8 - HI*HeI*k58 - HI*e*k1 - HI*e*k7 + HII*HM*k16 + \
             HII*e*k2 + HM*e*k14
    dHMdt = -H2II*HM*k19 - HI*HM*k8 + HI*e*k7 - HII*HM*k16 - HII*HM*k17 - HM*e*k14
    dHIIdt = -H2I*HII*k11 + H2II*HI*k10 - HI*HII*k9 + HI*HM*k8 + HI*HeI*k58 + \
             HI*e*k1 + 2*HI*k57 - HII*HM*k16 - HII*HM*k17 - HII*e*k2
    dHeIdt = -HeI*e*k3 + HeII*e*k4
    dHeIIdt = HeI*e*k3 - HeII*e*k4 - HeII*e*k5 + HeIII*e*k6
    dHeIIIdt = HeII*e*k5 - HeIII*e*k6
    dH2Idt = -H2I*HI*k13 - H2I*HII*k11 - H2I*e*k12 + H2II*HI*k10 + H2II*HM*k19 + HI*HM*k8 + 3*HI*k22
    dH2IIdt = H2I*HII*k11 - H2II*HI*k10 - H2II*HM*k19 - H2II*e*k18 + HI*HII*k9 + HII*HM*k17

    #"...since we have charge conservation, you don't *need* to solve the RHS for electron.
    #You can supply the RHS for e- as 0, and then recompute what the n_e is
    #every timestep instead (and use that as input to your RHS calculators)..."
    dedt = 0

    #dTdt = (T - state[9])
    dTdt = 0

    return np.array([dHIdt, dHMdt, dHIIdt, dHeIdt, dHeIIdt, dHeIIIdt,
                    dH2Idt, dH2IIdt, dedt, dTdt])




#Get user input if -v flag is specified
if(len(sys.argv) > 1):
    if(sys.argv[1] == "-v"):
        n_total = input("Input N where 10^N is the initial number of H, He, and H2...")
        h_ionized_frac = input("Input initial H ionization fraction...")
        he_ionized_frac = input("Input initial He ionization fraction...")
        h_mol_ionized_frac = input("Input initial molecular H ionization fraction...")
        T = input("Input temperature...")
        final_t = input("Input time to evolve to...")
        safety_factor = input("Input safety factor... ")
    else:
        print("Invalid flag!\nExiting...")
        exit(1)
#Otherwise default to some values
else:
    n_total = 5
    h_ionized_frac = -6
    he_ionized_frac = -5
    h_mol_ionized_frac = -2
    T = 150000
    final_t = 1000000
    safety_factor = 100000

#Initialize our state and calculate energy
n_HI_initial = 10**n_total * (1.0 - 10**h_ionized_frac)
n_HM_initial = 0
n_HII_initial = 10**n_total * 10**h_ionized_frac
n_HeI_initial = 10**n_total * (1.0 - 10**he_ionized_frac)
n_HeII_initial = 10**n_total * 10**he_ionized_frac
n_HeIII_initial = 0
n_H2I_initial = 10**n_total * (1.0 - 10**h_mol_ionized_frac)
n_H2II_initial = 10**n_total * 10**h_mol_ionized_frac
# n_e = n_HII + n_HeII + 2*n_HeIII + n_H2II - n_HM
n_e_initial = n_HII_initial + n_HeII_initial + n_H2II_initial

state_vector = np.array([n_HI_initial, n_HM_initial, n_HII_initial, \
                         n_HeI_initial, n_HeII_initial, n_HeIII_initial, \
                         n_H2I_initial, n_H2II_initial, n_e_initial, T])
energy = rr.energy_from_temp(state_vector, T)




#Setup and run integrator until final_t
integrator = sint.ode(updater)
integrator.set_initial_value(state_vector, t=0)
state_vector_values = []
ts = []
dt = final_t / safety_factor

ts.append(integrator.t)
state_vector_values.append(integrator.y)

while integrator.t < final_t:
    #Only need to compute these every time-step, not every call to updater()!
    #Cut my execution time in half :)
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, \
    k11, k12, k13, k14, k15, k16, k17, k18,  \
    k21, k22, k19, k57, k58 = calculate_reactionrates(T)

    integrator.integrate(integrator.t + dt)
    ts.append(integrator.t)

    #Recompute number of electrons every timestep
    #e = HII + HeII + 2*HeIII + H2II - HM
    integrator.y[8] = integrator.y[2] + integrator.y[4] + 2*integrator.y[5] \
                    + integrator.y[7] - integrator.y[1]

    #Recompute temperature from new number density
    T = rr.temperature(integrator.y, T, energy)
    integrator.y[9] = T

    state_vector_values.append(integrator.y)




#Graph state as a function of time
state_vector_values = np.array(state_vector_values)
ts = np.array(ts)
plt.loglog(ts, state_vector_values[:,0], label='HI')
plt.loglog(ts, state_vector_values[:,1], label='HM', color = 'b')
plt.loglog(ts, state_vector_values[:,2], label='HII')
plt.loglog(ts, state_vector_values[:,3], label='HeI')
plt.loglog(ts, state_vector_values[:,4], label='HeII')
plt.loglog(ts, state_vector_values[:,5], label='HeIII')
plt.loglog(ts, state_vector_values[:,6], label='H2I')
plt.loglog(ts, state_vector_values[:,7], label='H2II')
plt.loglog(ts, state_vector_values[:,8], label='e')
plt.loglog(ts, state_vector_values[:,9], label='T')

mpl.style.use('seaborn')

plt.xlabel("Time [s]")
plt.ylabel("N")
plt.legend()
plt.savefig("simulation.png")
