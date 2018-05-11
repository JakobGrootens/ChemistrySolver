import sympy

#Define our species symbolically
HI, HM, HII, HeI, HeII, HeIII, H2I, H2II, e = sympy.sympify(
"HI, HM, HII, HeI, HeII, HeIII, H2I, H2II, e")

#Define our reaction rates symbolically
k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, \
k11, k12, k13, k14, k15, k16, k17, k18,  \
k21, k22, k19, k57, k58 = sympy.sympify(
"k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, \
k11, k12, k13, k14, k15, k16, k17, k18,  \
k21, k22, k19, k57, k58")

#Define our reactions - sourced from Grackle paper
#  ---1:--       HI    + e   -> HII   + 2e
r1 = (HI + e), (HII + e + e), k1

#  ---2:--       HII   + e   -> H     + p
r2 = (HII + e), (HI), k2

#  ---3:--       HeI   + e   -> HeII  + 2e
r3 = (HeI + e), (HeII + e + e), k3

#  ---4:--       HeII  + e   -> HeI   + p
r4 = (HeII + e), (HeI), k4

#  ---5:--       HeII  + e   -> HeIII + 2e
r5 = (HeII + e), (HeIII + e + e), k5

#  ---6:--       HeIII + e   -> HeII  + p
r6 = (HeIII + e), (HeII), k6

#  ---7:--       HI    + e   -> HM    + p
r7 = (HI + e), (HM), k7

#  ---8:--       HM    + HI  -> H2I*  + e
r8 = (HM + HI), (H2I + HII), k8

#  ---9:--       HI    + HII -> H2II  + p
r9 = (HI + HII), (H2II), k9

#  ---10--       H2II  + HI  -> H2I*  + HII
r10 = (H2II + HI), (H2I + HII), k10

#  ---11--       H2I   + HII -> H2II  + H
r11 = (H2I + HII), (H2II + HI), k11

#  ---12--       H2I   + e   -> 2HI   + e
r12 = (H2I + e), (HI + HI + e), k12

#  ---13--       H2I   + H   -> 3H
r13 = (H2I + HI), (HI + HI + HI ), k13

#  ---14--       HM    + e   -> HI    + 2e
r14 = (HM + e), (HI + e + e), k14

#  ---15--       HM    + HI  -> 2H    + e
r15 = (HM + HI), (HI + HI + e), k15

#  ---16--       HM    + HII -> 2HI
r16 = (HM + HII), (HI + HI), k16

#  ---17--       HM    + HII -> H2II  + e
r17 = (HM + HII), (H2II + e), k17

#  ---18--       H2II  + e   -> 2HI
r18 = (H2II + e), (HI + HI), k18

#  ---19--       H2II  + HM  -> HI    + H2I
r19 = (H2II + HM), (HI + H2I), k19

#  ---21--       2H    + H2  -> H2I   + H2I
r21 = (HI + HI + H2I), (H2I + H2I), k21

#  ---22--       2H    + H   -> H2I   + H
r22 = (HI + HI + HI), (H2I + HI), k22

# ------- 57:   HI    + HI  -> HII   + HI    + e
r57 = (HI + HI), (HII + e + HI), k57

# ------- 58:   HI    + HeI -> HII   + HeI   + e
r58 = (HI + HeI), (HII + e + HeI), k58

#TODO: Get density dependant H + H + H and H + H + H2

all_reactions = [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10,
                r11, r12, r13, r14, r16, r17, r18, r19, r21, r22, r57, r58]

all_species = [HI, HM, HII, HeI, HeII, HeIII, H2I, H2II, e]


#Look for equations where a given species is
#formed - i.e. on the right hand side
def find_formation(species):
    f_reactions = []
    for reaction in all_reactions:
        if species in reaction[1].atoms():
            f_reactions.append(reaction)
    return f_reactions


#Look for equations where a given species is
#destroyed - i.e. on the left hand side
def find_destruction(species):
    d_reactions = []
    for reaction in all_reactions:
        if species in reaction[0].atoms():
            d_reactions.append(reaction)
    return d_reactions


#Create a right hand side for a given species
def get_rhs(species):
    dSdt = 0
    for lhs, rhs, coeff in find_formation(species):
        term = coeff
        for atom in list(lhs.atoms()):
            term *= atom
        dSdt += term
    for lhs, rhs, coeff in find_destruction(species):
        term = -coeff
        for atom in list(lhs.atoms()):
            term *= atom
        dSdt += term
    return dSdt

for species in all_species:
   print("{} RHS: {}\n".format(species, get_rhs(species)))
