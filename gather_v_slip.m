taxis = []
vcenter = []
for itt = 10:10:10000
	cmd = ['load snapshot',num2str(itt),'.mat'];
	eval(cmd);
	taxis = [taxis,t];
	vcenter = [vcenter,Vf(60)];
end
