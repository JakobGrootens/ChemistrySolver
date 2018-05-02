import sympy

#define our species symbolically
HI, HM, HII, HeI, HeII, HeIII, H2I, H2II, e = sympy.sympify(
"HI, HM, HII, HeI, HeII, HeIII, H2I, H2II, e")

#define our reaction rates symbolically
k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
k11, k12, k13, k14, k15, k16, k17, k18 = sympy.sympify(
"k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
 k11, k12, k13, k14, k15, k17, k18, k19")

#define our reactions - sourced from Grackle paper
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

r1 = (HI + de), (HII + de + de), k1
r2 = (HII + de), (HI), k2
r3 = (HeI + de), (HeII + de + de), k3
r4 = (HeII + de), (HeI), k4
r5 = (HeII + de), (HeIII + de + de), k5
r6 = (HeIII + de), (HeII), k6
r7 = (HI + HI), (HII + de + HI), k7
r8 = (HI + HeI), (HII + de + HeI), k8

#Additional nine species reactions
