function [Vf,tau]=NRsearch_NEW(fo,Vo,cca,ccb,Seff,tau,tauo,psi,...
    FaultZ,FaultVFree);

% limit direct effect by limiting the slip velocity 
Vw = 10^10;
fact = 1+(Vo/Vw)*exp(-psi);
% Vlimit = 10^10;
% lnVlVo = log(Vlimit/Vo);
% flimit = (fo+cca*lnVlVo+ccb*psi)/fact;
% taulimit = flimit*Seff;
% Vfrpos = (FaultZ*FaultVFree + tauo - taulimit)/(0.5*FaultZ);
% if (tau > 0 && Vfrpos > Vlimit)
%     Vf = Vfrpos;
%     tau = taulimit;
%     f = flimit;
%     return % Bypass N-R search
% end
% Vfrneg = (FaultZ*FaultVFree + tauo + taulimit)/(0.5*FaultZ);
% if (tau < 0 && Vfrneg < -Vlimit)
%     Vf = Vfrneg;
%     tau = -taulimit;
%     f = -flimit;
%     return % Bypass N-R search
% end

% N-R search (point-by-point) for tau if |Vf| < Vlimit:
cA = cca*Seff;
eps = 0.001*cA;
k=0;
delta=inf;
while (abs(delta) > eps)
    fa = fact*tau/(Seff*cca);
    help = -(fo+ccb*psi)/cca;
    help1 = exp(help+fa);
    help2 = exp(help-fa);
    Vf = Vo*(help1 - help2);
    Vfprime = fact*(Vo/cA)*(help1 + help2);
    delta = (FaultZ*FaultVFree - FaultZ*Vf + tauo - tau)/(1 + FaultZ*Vfprime);
    tau = tau + delta;
    k = k+1;
    if (abs(delta) > 10^10 || k == 100)
        k
        Vf
        tau
        Vfrpos
        error('N-R search fails to converge');
    end
end
% compute the updated Vf before exiting NR search
fa = fact*tau/(Seff*cca);
help = -(fo+ccb*psi)/cca;
help1 = exp(help+fa);
help2 = exp(help-fa);
Vf = Vo*(help1 - help2);


