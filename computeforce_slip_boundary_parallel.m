function [a] = computeforce_slip_boundary_parallel(iglob,W,Wl,H,Ht,d,slip,Iglob,coefint1, coefint2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    NGLOB = length(d);
    NGLL=5;
    [~,~,NEL] =size(iglob);
    a=zeros(NGLOB,2);
    ax = zeros(NGLOB,1);
    az = zeros(NGLOB,1);
    d(Iglob(:,1),2) = d(Iglob(:,2),2); % make a copy of displacement field.
    d(Iglob(:,1),1) = d(Iglob(:,2),1) + slip;
    local_fx =  zeros(NGLL,NGLL,NEL);
    local_fz =  zeros(NGLL,NGLL,NEL);
     dx = d(:,1);
     dz = d(:,2);
   parfor e=1:NEL              
        %switch to local (element) representation
        ig = iglob(:,:,e);

        local_x = dx(ig);
        local_z = dz(ig);
        %gradients wrt local variables (xi,eta)
    %    d_xi  = Ht*local;
    %    d_eta = local*H;
        %element contribution to internal forces
        %local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
        wloc = W(:,:,e);
        wlloc = Wl(:,:,e);
        
        dxxxx = H * ( (wloc + 2*wlloc) .* (Ht*local_x)) * coefint1; 
        dzzzz = ((wloc + 2*wlloc) .* (local_z*H)) * Ht * coefint2;
        dxxzz = H * (wlloc .* (local_z*H)); 
        dzzxx = (wlloc .* (Ht*local_x)) * Ht;
        
        dxzxz = (wloc .* (local_x*H)) * Ht * coefint2;
        dxzzx = (wloc .* (Ht*local_z)) * Ht;
        dzxxz = H * (wloc .* (local_x*H));
        dzxzx = H * (wloc .* (Ht*local_z)) * coefint1;
        
        
        local_fx(:,:,e) = dxxxx + dxxzz + dxzxz + dxzzx;
        local_fz(:,:,e) = dzzzz + dzzxx + dzxzx + dzxxz;
   end
   for e = 1:NEL
       ig = iglob(:,:,e);
        %assemble into global vector
        ax(ig) = ax(ig) + local_fx(:,:,e);
        az(ig) = az(ig) + local_fz(:,:,e);
    end
    a(:,1) = ax;
    a(:,2) = az;
    a(Iglob(:,2),2) = a(Iglob(:,1),2) +  a(Iglob(:,2),2); % combine vertical forces
    a(Iglob(:,2),1) = a(Iglob(:,1),1) +  a(Iglob(:,2),1); % combine horizontal forces
    a(Iglob(:,1),2) = 0;
    a(Iglob(:,1),1) = 0;
end

