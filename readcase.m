function [system, bus, branch, Dbus, OLTC, Dbranch, gen, avr, pss, gov, HVDC, fault, Y]=readcase(casename)
%% filename is casename
%% read case information and form the Y, B, and BB matrix
[bus, branch, Dbus, Dbranch, OLTC, gen, avr, pss, gov, DT, fault, system] = feval(casename);
%% build admittance matrices 变压器模型------K:1-----r+jx-------
jay = sqrt(-1);
nl = size(branch, 1);                          nb = size(bus, 1);
Zs = (branch(:,3) + jay * branch(:,4));  Zs(Zs==0)=1e-10;
Ys = 1 ./ Zs;
Bic = 1 .* branch(:,5);
Bjc = 1 .* branch(:,6);
tap = ones(nl, 1);
i = find(branch(:, 7));
tap(i) = branch(i, 7);
Ytt = Ys ;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;
Ytt = Ytt + jay*Bic;
Yff = Yff + jay*Bjc;
%% 生成交流系统的导纳矩阵
f = branch(:, 1);
t = branch(:, 2);
gl =find(f==t);
yg = sparse(size(branch,1),1);
yg(gl) = jay .* (branch(gl,5)+branch(gl,6));
Cf = sparse(f, 1:nl, ones(nl, 1), nb, nl);
Ct = sparse(t, 1:nl, ones(nl, 1), nb, nl);
Yg = Cf * spdiags(yg,0,nl,nl) * Cf';
Cf(:,gl)=0;
Ct(:,gl)=0;
pshunt = bus(:,11);
qshunt = bus(:,12);
Yq = sparse(1:nb,1:nb,pshunt+jay*qshunt,nb,nb);
Y = Cf * spdiags(Yff, 0, nl, nl) * Cf' + ...
	Cf * spdiags(Yft, 0, nl, nl) * Ct' + ...
	Ct * spdiags(Ytf, 0, nl, nl) * Cf' + ...
	Ct * spdiags(Ytt, 0, nl, nl) * Ct';
Y = Y + Yg+ Yq;
%% 发电机分区编号
gen(:,17) = bus(gen(:,2),15);
system.tcl = fault(6);
system.t0  = fault(5);

%% 直流系统初始化
if ~isempty(Dbus)
    HVDC.MVABASE = system.MVABASE;
    HVDC.UdcB   = Dbranch(:,9);
    HVDC.IdcB   = HVDC.MVABASE./HVDC.UdcB;
    HVDC.ZdcB   = HVDC.UdcB./HVDC.IdcB;
    if ~isempty(DT)
        HVDC.Tsr  = DT(Dbranch(:,1),2);
        HVDC.Kpr  = DT(Dbranch(:,1),3);
        HVDC.Kir  = DT(Dbranch(:,1),4);
        HVDC.Tsi  = DT(Dbranch(:,2),2);
        HVDC.Kpi  = DT(Dbranch(:,2),3);
        HVDC.Kii  = DT(Dbranch(:,2),4);
    else
        unit      = ones(size(Dbranch,1),1);
        HVDC.Kpr  = 0.04.*unit;
        HVDC.Kir  = 0.001.*unit;
        HVDC.Kpi  = 0.5.*unit;
        HVDC.Kii  = 5.0.*unit;
    end
else
    HVDC = [];
    return;
end
HVDC.bpaD = [];

%% 结束
return;

