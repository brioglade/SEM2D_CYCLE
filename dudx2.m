function [der] = dudx2(signal,H,Ht,wgll)
%DIFFUSE Summary of this function goes here
%   Detailed explanation goes here
u = H * diag(wgll) * Ht;
n = length(signal)/4;
signal(end+1) = signal(1);
f  = zeros(4*n+1,1);
for i = 1:1:n
    %i
    f(i*4-3:i*4+1) = f(i*4-3:i*4+1) + u*signal(i*4-3:i*4+1);
end
f(1) = f(1)+f(end);
m = zeros(4*n+1,1);
for i = 1:1:n
    m(i*4-3:i*4+1) = m(i*4-3:i*4+1) + wgll;
end
m(1) = m(1)+m(end);
der = f./m;
der = der(1:end-1);% second order derivative 
end

