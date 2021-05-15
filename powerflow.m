function  [bus, HVDC] = powerflow(system, bus, branch, Dbus, OLTC, Dbranch, HVDC, Y, gen, method)
if  method == 1
%     [bus, iter] = PQ_powerflow(system, bus, branch, Dbus, OLTC, Dbranch, Y);
    [bus, iter]=PQ_powerflow(bus,branch,Y);
else
    [bus, HVDC, iter] = NL_powerflow(system, bus, Dbus, OLTC, Dbranch, HVDC, Y);
end
%% 负荷模型
unit = ones(size(bus,1),1);
bus(:,19) = 1.0*unit;                bus(:,20) = 1.0*unit;
bus(:,21) = 0.0*unit;                bus(:,22) = 0.0*unit;
bus(:,23) = 0.0*unit;                bus(:,24) = 0.0*unit;
U       = bus(:,2);                  U2     = U.*U;
pload   = bus(:,6);                  qload  = bus(:,7);
p1      = bus(:,19);                 q1     = bus(:,20);
p2      = bus(:,21);                 q2     = bus(:,22);
p3      = bus(:,23);                 q3     = bus(:,24);
p0      = pload./(p1.*U2+p2.*U+p3);  q0     = qload./(q1.*U2+q2.*U+q3);
bus(:,24) = p0;                      bus(:,25) = q0;

%% 将非发电机的平衡节点的功率折算成负荷
Pe = bus(:,16);
Qe = bus(:,17);
vbus  = find(bus(:,4) > 0);
genbus = gen(:,2);
for i = 1:size(vbus,1)
    no = vbus(i);
    index = (find(genbus == no));
    if isempty(index)
        bus(no,4) = 0;
        bus(no,5) = 0;
        bus(no,6) = -Pe(no,1);
        bus(no,7) = -Qe(no,1);
    end
end

%% 换流器母线为PV节点
if ~isempty(HVDC)
    index = find(bus(HVDC.racbus,13) == 2);
    if ~isempty(index)
        Qac   = bus(HVDC.racbus(index),17);
        Qdc   = HVDC.Udr(index).*HVDC.Idr(index).*HVDC.tanfir(index);
        bus(HVDC.racbus(index),7) = - Qac - Qdc;
    end
    index = find(bus(HVDC.iacbus,13) == 2);
    if ~isempty(index)
        Qac   = bus(HVDC.iacbus(index),17);
        Qdc   = HVDC.Udi(index).*HVDC.Idi(index).*HVDC.tanfii(index);
        bus(HVDC.iacbus(index),7) = - Qac - Qdc;
    end
end

%% PV节点无发电机
% pvbus = find(abs(bus(:,4))>0 | abs(bus(:,5)>0) | abs(bus(:,8)>0) | abs(bus(:,9)>0));
% busno = gen(:,2);
% gvbus = find(abs(bus(busno,4))>0 | abs(bus(busno,5)>0));
% gvbus = gen(gvbus,2);
% index = setdiff(pvbus,gvbus);
% bus(index,6) = -Pe(index,1);
% bus(index,7) = -Qe(index,1);


% refbus = find(bus(:,13) == 3);
% Pe = bus(:,16);
% Qe = bus(:,17);
% index = find(bus(refbus,10) == 0);
% if ~isempty(index)
%     bus(refbus(index),6) = -Pe(refbus(index),1);
%     bus(refbus(index),7) = -Qe(refbus(index),1);
% end

%% 潮流输出
fprintf('Power Flow Result: iter = %2d\n',iter);
fprintf('\n  BusNum     Magn(p.u.)     Angle(°)        Pe(p.u.)        Qe(p.u.)');
fprintf('\n %5d %15.8f %15.8f %15.8f %15.8f',[bus(:,1), bus(:,2), bus(:,3)/pi*180, bus(:,16), bus(:,17)]');
fprintf('\n');
squeeze(iter);
%% 结束
return;