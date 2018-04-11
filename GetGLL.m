function [x,w,h]=GetGLL(ngll)
fid=fopen(sprintf('gll_%0.2u.tab',ngll));
data=fscanf(fid,'%f',[ngll,ngll+2]);
fclose(fid);
data
x=data(:,1)
w=data(:,2)
h=data(:,3:end)'
