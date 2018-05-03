import CalcReactionRates as rr
import scipy.integrate as sint
import matplotlib.pyplot as plt
import numpy as np


def updater(t, state):
    HI = state[0]; HM = state[1]; HII = state[2]; HeI = state[3]
    HeII = state[4]; HeIII = state[5]; H2I = state[6]; H2II = state[7]

    e = state[8]
    #Recompute number of electrons every timestep
    #e = HII + HeII + 2*HeIII + H2II - HM

    T = state[9]

    k1 = rr.k1(T); k2 = rr.k2(T); k3 = rr.k3(T)
    k4 = rr.k4(T); k5 = rr.k5(T); k6 = rr.k6(T)
    k7 = rr.k7(T); k8 = rr.k8(T); k9 = rr.k9(T)
    k10 = rr.k10(T); k11 = rr.k11(T); k12 = rr.k12(T)
    k13 = rr.k13(T); k14 = rr.k14(T); k15 = rr.k15(T)
    k16 = rr.k16(T); k17 = rr.k17(T); k18 = rr.k18(T)
    k19 = rr.k19(T);

    #Expressions made from RHSGenerator.py
    dHIdt =  H2I*HII*k11 + H2I*e*k12 - H2II*HI*k10 + H2II*HM*k19 + H2II*e*k18 - HI*HII*k9 - HI*HM*k8 - HI*e*k1 - HI*e*k7 + HII*HM*k16 + HII*e*k2 + HM*e*k14
    dHMdt = -H2II*HM*k19 - HI*HM*k8 + HI*e*k7 - HII*HM*k16 - HII*HM*k17 - HM*e*k14
    dHIIdt = -H2I*HII*k11 + H2II*HI*k10 - HI*HII*k9 + HI*HM*k8 + HI*e*k1 - HII*HM*k16 - HII*HM*k17 - HII*e*k2
    dHeIdt = -HeI*e*k3 + HeII*e*k4
    dHeIIdt = HeI*e*k3 - HeII*e*k4 - HeII*e*k5 + HeIII*e*k6
    dHeIIIdt = HeII*e*k5 - HeIII*e*k6
    dH2Idt = -H2I*HI*k13 - H2I*HII*k11 - H2I*e*k12 + H2II*HI*k10 + H2II*HM*k19 + HI*HM*k8
    dH2IIdt = H2I*HII*k11 - H2II*HI*k10 - H2II*HM*k19 - H2II*e*k18 + HI*HII*k9 + HII*HM*k17

    #"...since we have charge conservation, you don't *need* to solve the RHS for electron.
    #You can supply the RHS for e- as 0, and then recompute what the n_e is
    #every timestep instead (and use that as input to your RHS calculators)..."

    dedt = dHIIdt + dHeIIdt + 2*dHeIIIdt + dH2IIdt - dHMdt

    #TODO update temperature using gas laws
    dTdt = 0

    ret = np.array([dHIdt, dHMdt, dHIIdt, dHeIdt, dHeIIdt, dHeIIIdt,
                    dH2Idt, dH2IIdt, dedt, dTdt])

    #print(state)
    return ret


n_total = 5
e_frac = -2
T = np.log10(150000)
final_t = 10000
safety_factor = 1000

n_HI_initial = 10**n_total * (1.0 - 10**e_frac)
n_HM_initial = 0
n_HII_initial = 10**n_total * 10**e_frac
n_HeI_initial = 10**n_total *.5
n_HeII_initial = 0
n_HeIII_initial = 0
n_H2I_initial = 0
n_H2II_initial = 0
# n_e = n_HII + n_HeII + 2*n_HeIII + n_H2II - n_HM
n_e_initial = 10**n_total * 10**e_frac
state_vector = np.array([n_HI_initial, n_HM_initial, n_HII_initial, \
                         n_HeI_initial, n_HeII_initial, n_HeIII_initial, \
                         n_H2I_initial, n_H2II_initial, n_e_initial, 10**T])

integrator = sint.ode(updater)
integrator.set_initial_value(state_vector, t=0)
state_vector_values = []
ts = []
dt = final_t / safety_factor

ts.append(integrator.t)
state_vector_values.append(integrator.y)
while integrator.t < final_t:
    integrator.integrate(integrator.t + dt)
    ts.append(integrator.t)
    state_vector_values.append(integrator.y)


state_vector_values = np.array(state_vector_values)
ts = np.array(ts)
plt.loglog(ts, state_vector_values[:,0], label='HI')
plt.loglog(ts, state_vector_values[:,1], label='HM')
plt.loglog(ts, state_vector_values[:,2], label='HII')
plt.loglog(ts, state_vector_values[:,3], label='HeI')
plt.loglog(ts, state_vector_values[:,4], label='HeII')
plt.loglog(ts, state_vector_values[:,5], label='HeIII')
plt.loglog(ts, state_vector_values[:,6], label='H2I')
plt.loglog(ts, state_vector_values[:,7], label='H2II')
plt.loglog(ts, state_vector_values[:,8], label='e')


plt.xlabel("Time [s]")
plt.ylabel("n")
plt.legend()

plt.savefig("simulation.png")
