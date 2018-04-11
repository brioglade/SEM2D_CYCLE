% Quasi-static + dynamic modeling of earthquake nucleation using 2D anti-plane SEM
% Variable time stepping has added
% Jacobi preconditioner added
% Output variables saved in one single file

%%%%%%%%%%%%%%%%%%% modified smooth ICs %%%%%%%%%%%%%%%%%%%%%%%%%%%
% a = 0.0144
% b = 0.0191
% L = 42 mm

%------------------------------------------
clear all; close all;
bound = 50;
%%%%%%%%%%% Initial Conditions and state variable evolution %%%%%%%%%%%%%%%

% If IDinitCond = 1, use SCEC initial conditions
% If IDinitCond = 2, use smooth initial conditions
IDintitCond = 2;

% If IDstate = 1, compute psi(t+dt) = psi(t) + dt * dpsi(t)/dt
% If IDstate = 2, compute psi(t+dt) by integration with constant V
% If IDstate = 3, compute psi(t+dt) of slip law by integration with constant V
IDstate = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1: SPECTRAL ELEMENT MESH GENERATION

%**** Set here the parameters of the square box domain and mesh : ****
distN = 1000;
NSIZE = 2;
LX = NSIZE*45e3/distN;
LY = NSIZE*30e3/distN;

yr2sec = 365*24*60*60;
tevent = 197392769.833719193;  % target event time (s)

if IDintitCond == 1
    Total_time = 4.5;    % total time (sec) for SECE initial conditions
elseif IDintitCond == 2 || IDintitCond == 3
    Total_time = 100*yr2sec; %22*yr2sec;  %250*yr2sec;
    %Total_time = 39.5;  % total time (sec) for soothe IC
end
   
%%% Note: I use st_node_space = '3_16' in the paper
%st_node_space = '1_8';  NELX = 90; NELY = 60; P = 4;
%st_node_space = '3_16'; NELX = 60; NELY = 40; P = 4;
%st_node_space = '1_4';  NELX = 45; NELY = 30; P = 4; 
%st_node_space = '3_8';  NELX = 30; NELY = 20; P = 4; 
st_node_space = '3_4';  NELX = 15; NELY = 10; P = 4; 

NELX = NELX*NSIZE;
NELY = NELY*NSIZE;

dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element

%[iglob,x,y] = MeshBox(LX,LY,NELX,NELY,NGLL);

XGLL = GetGLL(NGLL);   % local cordinate of GLL quadrature points

iglob = zeros(NGLL,NGLL,NELX*NELY);	% local to global index mapping
nglob = (NELX*(NGLL-1)+1)*(NELY*(NGLL-1)+1);	% number of global nodes

x     = zeros(nglob,1);		% coordinates of GLL nodes
y     = zeros(nglob,1);	

e=0;
last_iglob = 0;
igL = reshape([1:NGLL*(NGLL-1)],NGLL-1,NGLL);
igB = reshape([1:NGLL*(NGLL-1)],NGLL,NGLL-1);
igLB = reshape([1:(NGLL-1)*(NGLL-1)],NGLL-1,NGLL-1);
xgll = repmat( 0.5*(1+XGLL) , 1,NGLL);
ygll = dye*xgll';
xgll = dxe*xgll;

for ey=1:NELY, 
for ex=1:NELX, 

  e = e+1;

 % Take care of redundant nodes at element edges :
  if e==1  % first element: bottom-left
    ig = reshape([1:NGLL*NGLL],NGLL,NGLL);    
  else
    if ey==1 	%  elements on first (bottom) row
      ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
      ig(2:NGLL,:) = last_iglob + igL; 		% the rest
    elseif ex==1 % elements on first (left) column
      ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
      ig(:,2:NGLL) = last_iglob + igB; 		% the rest
    else 	% other elements
      ig(1,:) = iglob(NGLL,:,e-1); 		% left edge
      ig(:,1) = iglob(:,NGLL,e-NELX); 		% bottom edge
      ig(2:NGLL,2:NGLL) = last_iglob + igLB;
    end
  end
  iglob(:,:,e) = ig;
  last_iglob = ig(NGLL,NGLL);

 % Global coordinates of the computational (GLL) nodes
  x(ig) = dxe*(ex-1)+xgll;
  y(ig) = dye*(ey-1)+ygll;

   
end
end

x = x-LX/2;
nglob = length(x);

RHO = 2670;
VS1 = 0.8*3464;  % modified for layer case
VS2 = 3464;
VP1 = 0.8*6000;
VP2 = 6000;
ETA = 0;
THICK = 30;  % Thickness of the first layer with VS1 (m)
            % Make sure THICK is integer multiple of dye 
%------------------------------------------
% STEP 2: INITIALIZATION

% 
[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

W     = zeros(NGLL,NGLL,NEL);	% for internal forces
Wl    = zeros(NGLL,NGLL,NEL);
M     = zeros(nglob,1);		    % global mass matrix, diagonal
rho   = zeros(NGLL,NGLL);	    % density will not be stored
mu    = zeros(NGLL,NGLL);	    % shear modulus will not be stored
muMax = 0;                    % Used for variable time stepping

%**** Set here the parameters of the time solver : ****
%NT = 347;                     % number of timesteps
CFL = 0.6;                     % stability number = CFL_1D / sqrt(2)
%********

dt = Inf;  			% timestep (set later)

% For this simple mesh the global-local coordinate map (x,y)->(xi,eta)
% is linear, its jacobian is constant
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;
coefint1 = jac/dx_dxi^2 ;
coefint2 = jac/dy_deta^2 ;

% FOR EACH ELEMENT ...
for ey=1:NELY,
    for ex=1:NELX,

        e = (ey-1)*NELX+ex;
        ig = iglob(:,:,e);
        
        %**** Set here the physical properties of the heterogeneous medium : ****
        if (ey*dye <= -10*THICK)  % modified for layer case
        rho(:,:) = RHO;
        mu(:,:)  = RHO* VS1^2;
        ld(:,:)  = RHO* VP1^2 - 2* mu(:,:);
        else
        rho(:,:) = RHO;
        mu(:,:)  = RHO* VS2^2;   
        ld(:,:)  = RHO* VP2^2 - 2* mu(:,:);
        end
        if muMax < max(max(mu)); muMax = max(max(mu)); end;
        
        %********
        
        % Diagonal mass matrix
        M(ig) = M(ig) + wgll2 .*rho *jac;
        
        % Local contributions to the stiffness matrix K
        %  WX(:,:,e) = wgll2 .* mu *jac/dx_dxi^2;
        %  WY(:,:,e) = wgll2 .* mu *jac/dy_deta^2;
        W(:,:,e) = wgll2 .* mu;
        Wl(:,:,e) = wgll2 .* ld;
        % The timestep dt is set by the stability condition
        % dt = CFL*min(dx/vs)
        vs = sqrt(mu./rho);
        if dxe<dye
            vs = max( vs(1:NGLL-1,:), vs(2:NGLL,:) );
            dx = repmat( diff(xgll)*0.5*dxe ,1,NGLL);
        else
            vs = max( vs(:,1:NGLL-1), vs(:,2:NGLL) );
            dx = repmat( diff(xgll)'*0.5*dye ,NGLL,1);
        end
        dtloc = dx./vs;
        dt = min( [dt dtloc(1:end)] );

    end
end  %... of element loop
dt = CFL*dt;
dtmin = dt;

tmax = Total_time;
dtmax = 50*24*60*60/distN*100; % 5 days

if ETA, dt=dt/sqrt(1+2*ETA); end  % dt modified slightly for damping
NT = ceil(Total_time/dtmin)       % estimate the max possible time step
half_dt = 0.5*dtmin;
half_dt_sq = 0.5*dtmin^2;

%-- Initialize kinematic fields, stored in global arrays
d = zeros(nglob,2);  %xy plane d(:,1) = ux, d(:,2) = uy;
v = zeros(nglob,2);
v(:,1) = 1/2*10^-3;
vPre = zeros(nglob,2);
a = zeros(nglob,2);

Vpl = 2*10^-3/yr2sec;    % imposed loading 2 mm/yr 

%-- Absorbing boundaries (first order):
% Left
ng = NELY*(NGLL-1)+1;
BcLeftIglob = zeros(ng,1);
BcLeftC = zeros(ng,1);
for ey=1:NELY,
    ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
    e = (ey-1)*NELX+1;
    BcLeftIglob(ip) = iglob(1,1:NGLL,e);
    if (ey*dye <= THICK)  % modified for layer case
        impedance = RHO*VS1;
    else
        impedance = RHO*VS2;         
    end
    BcLeftC(ip) = BcLeftC(ip) + dy_deta*wgll*impedance ;
end
% Right
ng = NELY*(NGLL-1)+1;
BcRightIglob = zeros(ng,1);
BcRightC = zeros(ng,1);
for ey=1:NELY,
    ip = (NGLL-1)*(ey-1)+[1:NGLL] ;
    e = (ey-1)*NELX+NELX;
    BcRightIglob(ip) = iglob(NGLL,1:NGLL,e);
    if (ey*dye <= THICK)  % modified for layer case
        impedance = RHO*VS1;
    else
        impedance = RHO*VS2;         
    end
    BcRightC(ip) = BcRightC(ip) + dy_deta*wgll*impedance ;
end
% Top
ng = NELX*(NGLL-1)+1;
BcTopIglob = zeros(ng,1);
BcTopC = zeros(ng,1);
for ex=1:NELX,
    ip = (NGLL-1)*(ex-1)+[1:NGLL] ;
    e = (NELY-1)*NELX+ex;
    BcTopIglob(ip) = iglob(1:NGLL,NGLL,e);
%    if (ex < NELX/2)          % modified for layer case   
%    impedance = RHO*VS1;
%    else
    impedance = RHO*VS2;
%    end
    BcTopC(ip) = BcTopC(ip) + dx_dxi*wgll*impedance ;
end

Mq = M;
% % The mass matrix needs to be modified at the boundary
% % for the IMPLICIT treatment of the term C*v.
% Fortunately C is diagonal.
M(BcLeftIglob)  = M(BcLeftIglob)  +half_dt*BcLeftC;
M(BcRightIglob) = M(BcRightIglob) +half_dt*BcRightC;
M(BcTopIglob) = M(BcTopIglob) + half_dt*BcTopC;


%-- DYNAMIC FAULT at bottom boundary
FaultNglob = NELX*(NGLL-1)+1;
FaultIglob = zeros(FaultNglob,1);
FaultB = zeros(FaultNglob,1);
for ex = 1:NELX,
    ip = (NGLL-1)*(ex-1)+[1:NGLL];
    e = ex;
    FaultIglob(ip) = iglob(1:NGLL,1,e);
    FaultB(ip) = FaultB(ip) + dx_dxi*wgll;    
end
%FaultZ = M(FaultIglob)./FaultB /half_dt;
FaultZ = M(FaultIglob)./FaultB /half_dt * 0.5;  % times 0.5 due to the symmetry
FaultX = x(FaultIglob);

Seff = repmat(120*10^6,FaultNglob,1);  % effective normal stress
fo = repmat(0.6,FaultNglob,1);         % reference friction
cca = repmat(0.0144,FaultNglob,1);      % constitutive parameter a
ccb = repmat(0.0191,FaultNglob,1);      % constitutive parameter b
Vo = repmat(10^-6,FaultNglob,1);       % reference velocity Vo
xLf = repmat(2*0.042/distN,FaultNglob,1);      % L (D_c) = 84 um
FaultC = repmat(0,FaultNglob,1);       % used only for integrated state variable 
Vf1 = repmat(0,FaultNglob,1);
Vf2 = repmat(0,FaultNglob,1);
Vf  = repmat(0,FaultNglob,1);
psi1 = repmat(0,FaultNglob,1);
psi2 = repmat(0,FaultNglob,1);
tau1 = repmat(0,FaultNglob,1);
tau2 = repmat(0,FaultNglob,1);
tau3 = repmat(0,FaultNglob,1);
tauNR = repmat(0,FaultNglob,1);
tauAB = repmat(0,FaultNglob,1);         % USED FOR QUASISTATIC

if IDintitCond == 1
     
      
elseif IDintitCond == 2
    
    %-- Initial conditions smooth in time and space
    tauoBack = 70*10^6;
    tauo = repmat(tauoBack,FaultNglob,1);
    width = 2*3e3/distN;
    isel = find(abs(FaultX)<=width/2);
    Peaktauo = 81.6*10^6;
    Amplitude = (Peaktauo-tauoBack)/2;
    tauo(isel) = (Peaktauo+tauoBack)/2 ...
        + Amplitude*cos(2*pi*FaultX(isel)/width);
    isel2 = find(abs(FaultX)>10e3/distN);
    
    ccbOut = 0.0097;
    ccbIn = 0.0191; 
    
    ccb(isel2) = ccbOut;
    Amplitude2 = (ccbIn + ccbOut)/2;
    Amplitude3 = (ccbIn - ccbOut)/2;
    isel3 = find(FaultX>=8e3/distN&FaultX<=(8e3/distN+width/2));
    ccb(isel3)=Amplitude2 + Amplitude3*cos(2*pi*(FaultX(isel3)-8e3)/width);
    isel4 = find(FaultX<=-8e3/distN&FaultX>=-(8e3/distN+width/2));
    ccb(isel4)=Amplitude2 + Amplitude3*cos(2*pi*(FaultX(isel4)+8e3)/width);
    tau = repmat(0,FaultNglob,1);
    psi = tauo./(Seff.*ccb) - fo./ccb - (cca./ccb).*log(2*v(FaultIglob)./Vo);
    psi0 = psi;

end

%d = d_store;
%v = v_store;
%psi = psi_store;
%psi0 = psi;

if ETA,  % Kelvin-Voigt viscosity
  %NEL_ETA = min( NELX, ceil(L_BARRIER/dxe)+2 );
  NEL_ETA = NELX;
  x1 = 0.5*(1+xgll');
  eta_taper = exp(-pi*x1.^2); 
  eta = ETA*dt *repmat(eta_taper, NGLL,1 );
else
  NEL_ETA = 0;
end

%-- initialize data for output seismograms
%**** Set here receiver locations : ****
OUTxseis = [-15:3:0]';     		    % x coord of receivers
OUTnseis = length(OUTxseis);		% total number of receivers
OUTyseis = repmat(15,OUTnseis,1);	% y coord of receivers
%********
%receivers are relocated to the nearest node
%OUTdseis = distance between requested and relocated receivers
[OUTxseis,OUTyseis,OUTiglob,OUTdseis] = FindNearestNode(OUTxseis,OUTyseis,x,y);
kkseis=1;
%OUTv = zeros(OUTnseis,NT);

%-- initialize data for output snapshots
OUTindx = Init2dSnapshot(iglob);

% time is a array stored time at each time step
t = 0;
it = 0;
dtincf = 1.2;
dtpre = dt;
gamma = pi/4;
% average node spacing
hcell = LX/(FaultNglob-1);
Ximax = 0.5;   % value from NL00 (15ab)
Xithf = 1;
trec = 0;
Vthres = 0.01;  % 1 cm/s
slipstart = 0;
ievb = 0;
ieva = 0;
ntvsx=0;
tvsx = 0.5*yr2sec;
tvsxinc = tvsx;
nevne = 0;
tevneinc = 0.001;
Vevne = Vthres;

% compute XiLf for each fault node
for j=1:FaultNglob
    % Compute time-restricting parameters as in section 4, NL00
    expr1=-(cca(j) - ccb(j))/cca(j);
    expr2=gamma*muMax/hcell*xLf(j)/(cca(j)*Seff(j));
    ro=expr2-expr1;
    if (0.25*ro*ro-expr2) >= 0 
        Xith(j)=1/ro;
    else    
        Xith(j)=1-expr1/expr2;
    end 
    % For each node, compute the slip that the node cannot
    % exceed in one time step; store in array XiLf(FaultNglob)
    if (Xithf*Xith(j) > Ximax) 
        XiLf(j)=Ximax*xLf(j);
    else
        XiLf(j)=Xithf*Xith(j)*xLf(j);
    end    
end

% time-related variables added
OUTt = 0.5;
q = 1;
OUTtGo = 0;
OUTtCount = 1;

% OUTPUT field quantities at several locations on fault
OutLoc1 = 0e3/distN;                % 0 km point
[OUTxLoc1,OUTyLoc1,OUTiglobLoc1] = FindNearestNode(OutLoc1,0,x,y);
FaultLoc1 = round((OutLoc1+LX/2)*FaultNglob/LX);

OutLoc2 = 3e3/distN;                % 3 km point
[OUTxLoc2,OUTyLoc2,OUTiglobLoc2] = FindNearestNode(OutLoc2,0,x,y);
FaultLoc2 = round((OutLoc2+LX/2)*FaultNglob/LX);

OutLoc3 = 6e3/distN;                % 6 km point
[OUTxLoc3,OUTyLoc3,OUTiglobLoc3] = FindNearestNode(OutLoc3,0,x,y);
FaultLoc3 = round((OutLoc3+LX/2)*FaultNglob/LX);

OutLoc4 = 9e3/distN;                % 9 km point
[OUTxLoc4,OUTyLoc4,OUTiglobLoc4] = FindNearestNode(OutLoc4,0,x,y);
FaultLoc4 = round((OutLoc4+LX/2)*FaultNglob/LX);

% center of model domain              
[OUTxLoc5,OUTyLoc5,OUTiglobLoc5] = FindNearestNode(0,30,x,y);

disp('Total number of nodes on fault = ');
disp(num2str(FaultNglob));
disp('Average node spacing = ');
disp(LX/1000/(FaultNglob-1));

fprintf('dt = %1.17f \n',dt);
pp = 1;

jj = 1;
jjj = 1;
for ii = 1:nglob
    if ii == FaultIglob(jj)
        jj = jj + 1; 
        if jj > length(FaultIglob)
            jj = length(FaultIglob);
        end
    else
       FaultNIglob(jjj,1) = ii;  % find nodes that do not belong to fault
       jjj = jjj + 1; 
    end
end

r = zeros(nglob,2);	
beta = zeros(nglob,1);	
alpha = zeros(nglob,1);	
p = zeros(nglob,2);	
F = zeros(nglob,2);	
dPre = zeros(nglob,2);
vPre = zeros(nglob,2);
dd = zeros(nglob,2);
dacum = zeros(nglob,2);

dnew = zeros(length(FaultNIglob),2);	
Fnew = zeros(length(FaultNIglob),2);	
anew = zeros(length(FaultNIglob),2);	
rnew = zeros(length(FaultNIglob),2);	
pnew = zeros(length(FaultNIglob),2);	 

%------------------------------------------
% STEP 3: SOLVER  M*a = -K*d +F
% Explicit Newmark scheme with
% alpha=1, beta=0, gamma=1/2
%
isolver = 1;   % initially, 1 for quasistatic, 2 for dynamic
go_snap = 0;
go_snapDY = 0;

% compute the diagonal of K
Kdiagx = zeros(nglob,1);	
Kdiagz = zeros(nglob,1);
Kdiag = zeros(nglob,2);

Klocdiagx = zeros(NGLL,NGLL);
Klocdiagz = zeros(NGLL,NGLL);

for e=1:NEL, % FOR EACH ELEMENT ...
    ig = iglob(:,:,e);
    wloc = W(:,:,e);
    wlloc = Wl(:,:,e);
    Klocdiagx(:,:) = 0;
    Klocdiagz(:,:) = 0;

    for k = 1:NGLL
        for j = 1:NGLL
            Klocdiagx(k,j)=Klocdiagx(k,j)+sum(coefint1*H(k,:)'.*(wlloc(:,j)+2*wloc(:,j)).*Ht(:,k)...
                    +coefint2*(wloc(k,:).*H(j,:))'.*Ht(:,j));
            Klocdiagz(k,j)=Klocdiagz(k,j)+sum(coefint1*H(k,:)'.*(wlloc(:,j)).*Ht(:,k)...
                    +coefint2*((2*wloc(k,:)+wlloc(k,:)).*H(j,:))'.*Ht(:,j));
        end
    end
    Kdiagx(ig) = Kdiagx(ig) + Klocdiagx(:,:);
    Kdiagz(ig) = Kdiagz(ig) + Klocdiagz(:,:);
end

Kdiag(:,1) = Kdiagx;
Kdiag(:,2) = Kdiagz;
diagKnew = Kdiag(FaultNIglob,:);

v(:,1) = v(:,1) - 0.5*Vpl;
Vf = 2*v(FaultIglob,:);
iFBC = find(abs(FaultX)>=bound*10^3/distN);
NFBC = length(iFBC);
Vf(iFBC,:) = 0;

jj = 1;
FaultIglobBC = zeros(NFBC,1);
for ex = 1:NELX
    for k = 1:NGLL
       if abs(x(iglob(k,1,ex))) >= bound*10^3/distN 
          FaultIglobBC(jj) = iglob(k,1,ex);
          jj = jj + 1;
       end
    end        
end
v(FaultIglobBC,:) = 0;


testX = 0*x;
testZ = 0*testX;
testX(FaultIglob) = cos(0.2*x(FaultIglob));
F = [testZ,testX];
dnew = zeros(length(FaultNIglob),2);
[dnew,~]=myPCGnew3(coefint1,coefint2,diagKnew,dnew,F,FaultIglob,...
        FaultNIglob,H,Ht,iglob,NEL,nglob,W,Wl,x,y);
    
    
    
    
    
    
    