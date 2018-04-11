 function [dnew,n] = myPCGnew2(coefint1,coefint2,diagKnew,dnew,F,Iglob,NIglob,H,Ht,iglob,NEL,nglob,W);

% input: 

% output: dnew, updated displacement
%         n,    number of PCG iteration

a(nglob,:) = 0;
dd(nglob,:) = 0;
p(nglob,:) = 0;

% first compute F = K21 d1, where d1 is incremental disp at fault node
for e=1:NEL,
    %switch to local (element) representation
    ig = iglob(:,:,e);
    local = F(ig);
    %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;
    d_eta = local*H;
    %element contribution to internal forces
    local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
    %assemble into global vector
    a(ig) = a(ig) + local;
end                 
Fnew = -a(NIglob);      % make it F = -K21 d1 
    
dd(NIglob) = dnew;
dd(Iglob) = 0;
% Solve d = K22^-1 * F by conjugate gradient
% first compute residual based on initial guess, dd
a(:) = 0;
for e=1:NEL,
    %switch to local (element) representation
    ig = iglob(:,:,e);
    local = dd(ig);
    %gradients wrt local variables (xi,eta)
    d_xi  = Ht*local;
    d_eta = local*H;
    %element contribution to internal forces
    local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
    %assemble into global vector
    a(ig) = a(ig) + local;
end
anew = a(NIglob);

rnew = Fnew - anew;                   % initial residual F - K22 d2
znew = rnew./diagKnew;
pnew = znew;
p(:) = 0;
p(NIglob) = pnew;
   
for n = 1:4000
    anew(:)=0;
    a(:) = 0;
    for e=1:NEL,
        %switch to local (element) representation
        ig = iglob(:,:,e);
        local = p(ig);
        %gradients wrt local variables (xi,eta)
        d_xi  = Ht*local;
        d_eta = local*H;
        %element contribution to internal forces
        local = coefint1*H*( W(:,:,e).*d_xi ) + coefint2*( W(:,:,e).*d_eta )*Ht ;
        %assemble into global vector
        a(ig) = a(ig) + local;
    end
    anew = a(NIglob);
    alpha = transpose(znew)*rnew/(transpose(pnew)*anew);
    dnew = dnew + alpha*pnew;
    r_old = rnew;
    z_old = znew;
    rnew = r_old - alpha*anew;
    znew = rnew./diagKnew;
    beta = transpose(znew)*rnew/(transpose(z_old)*r_old);
    pnew = znew + beta*pnew;
    p(:) = 0;
    p(NIglob) = pnew;
    %ResStore(n) = norm(rnew)/norm(Fnew);
    if norm(rnew)/norm(Fnew) < 10^-5
        %'converged'
        %figure;
        %loglog([1:n],ResStore(1:n),'b.-');
        %axis([10^0 10^3 10^-4 10^0]);
        break
    end
    if n == 4000 || norm(rnew)/norm(Fnew) > 10^10
        'CG not converge'
%         figure;
%         loglog([1:n],ResStore(1:n),'b.-');
%         axis([10^0 10^3 10^-4 10^0]);
        return
    end
end
