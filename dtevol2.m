function [dt,trans]=dtevol2(dt,dtmax,dtmin,dtincf,XiLf,FaultNglob,NFBC,Vf,isolver);

% input: arrays - Vf, XiLf
%        real - dt, dtmax, dtmin, dtincf
%        integer - FaultNglob
%
%
%
% output: dt

if isolver == 1 % for quasi-static solver

    dtnx = dtmax;          % dinitial value of dt

    % The next loop adjusts the time step according to cell velocities
    % and allowed slip in one time step:
    for iF=1:FaultNglob-NFBC;
        i = iF + NFBC/2;
        absV = abs(Vf(i,1));
        if (absV*dtmax > XiLf(i))   % velocity is not too small
            dtcell = XiLf(i)/(absV);  % then see whether cell i restrics the time step
            if (dtcell < dtnx)      % the factor 40 is added by ykaneko
                dtnx = dtcell;
            end
        end
    end

    if (dtmin > dtnx)      % dtmin is minmum time step allowed
        dtnx = dtmin;
    end

    if (dtnx > dtincf*dt)  % limit dt increase
        dtnx = dtincf*dt;
    end
   
    dt = dtnx;             % choose time step to be dt;
    
    if(dt == dtmin)
        trans = true;
    else
        trans = false;
    end
elseif isolver == 2        % for dynamic solver

    dt = dtmin;
    trans = false;

end
