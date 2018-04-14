 function [dnew,n] = myPCGnew5(coefint1,coefint2,diagKnew,dnew,F,Iglob,NIglob,H,Ht,iglob,NEL,nglob,W,Wl,FixBoundary,x,y);

% input: 

% output: dnew, updated displacement
%         n,    number of PCG iteration

a(nglob,2) = 0;
dd(nglob,2) = 0;
p(nglob,2) = 0;

% first compute F = K21 d1, where d1 is incremental disp at fault node
a = computeforce_combine_force(iglob,W,Wl,H,Ht,F,Iglob,coefint1,coefint2);               
a(FixBoundary,:) = 0;
a(Iglob,1) = 0;
Fnew = -a;      % make it F = -K21 d1 

dd(NIglob,:) = dnew(NIglob,:);
dd(Iglob,2) = dnew(Iglob,2);
dd(Iglob,1) = 0;
% Solve d = K22^-1 * F by conjugate gradient
% first compute residual based on initial guess, dd
a = computeforce_combine_force(iglob,W,Wl,H,Ht,dd,Iglob,coefint1,coefint2);% combine vertical stresses 
a(FixBoundary,:) = 0;
anew(NIglob,:) = a(NIglob,:);    %For fixed boundary
anew(Iglob,2) = a(Iglob,2);
anew(Iglob,1) = 0;

rnew = Fnew - anew;                   % initial residual F - K22 d2
znew = rnew./diagKnew;
pnew = znew;
p(:) = 0;
p(NIglob,:) = pnew(NIglob,:);
p(Iglob,2) = pnew(Iglob,2);
   
for n = 1:4000
    a = computeforce_combine_force(iglob,W,Wl,H,Ht,p,Iglob, coefint1,coefint2);
    a(FixBoundary,:) = 0;
    anew(NIglob,:) = a(NIglob,:);
    anew(Iglob,2) = a(Iglob,2);
    
    alpha = sum(diag(transpose(znew)*rnew))/sum(diag((transpose(pnew)*anew)));
    dnew = dnew + alpha*pnew;
    r_old = rnew;
    z_old = znew;
    rnew = r_old - alpha*anew;
    znew = rnew./diagKnew;
    beta = sum(diag(transpose(znew)*rnew))/sum(diag(transpose(z_old)*r_old));
    pnew = znew + beta*pnew;
    p(:,:) = 0;
    p(NIglob,:) = pnew(NIglob,:);
    p(Iglob,2) = pnew(Iglob,2);
    %ResStore(n) = norm(rnew)/norm(Fnew);
    if norm(rnew)/norm(Fnew) < 10^-6 && norm(rnew) < 1.0
        %'converged'
        %figure;
        %loglog([1:n],ResStore(1:n),'b.-');
        %axis([10^0 10^3 10^-4 10^0]);
        
        display(['norm:',num2str(norm(rnew)),',iterations:',num2str(n)]);
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
