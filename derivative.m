syms pH_2 pH_p pH(pH_2,pH_p) k1 k2(pH_2,pH_p) k3(pH_2,pH_p) nH(pH_2,pH_p) nH_2(pH_2)

jacobian(k1*pH(pH_2,pH_p)*nH(pH_2,pH_p) - k2(pH_2,pH_p)*pH_2*nH(pH_2,pH_p) - k3(pH_2,pH_p)*pH_2*nH_2(pH_2), [pH_2, pH_p])