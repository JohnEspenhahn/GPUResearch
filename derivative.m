% syms pH_p t N xC xSi M_h ne(pH_p)
% psi = symfun((1.13*sqrt(t))/ne(pH_p), pH_p);
% diff(N*(xC + xSi + pH_p / (N*1.36*M_h)), pH_p)

syms pH2 pH_p pH(pH2,pH_p) k1 k2(pH2,pH_p) k3(pH2,pH_p) nH(pH2,pH_p) nH_2(pH2) ne(pH_p) k6 k7 k8(pH_p)
% a = jacobian(k1*pH(pH2,pH_p)*nH(pH2,pH_p) - k2(pH2,pH_p)*pH2*nH(pH2,pH_p) - k3(pH2,pH_p)*pH2*nH_2(pH2), [pH2, pH_p])
b = jacobian(k6*pH(pH2,pH_p)*ne(pH_p) - k7*pH_p*ne(pH_p) - k8(pH_p)*pH_p*ne(pH_p), [pH2, pH_p])
% diff((1e-14*12.25)/(1 + 8.074e-6*psi^1.378*(1 + 5.087e2*t^1.586e-2*psi^(-0.4723 - 1.102e-5*log(t)))), pH_p)