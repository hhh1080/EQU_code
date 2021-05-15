clear;
clc;
close all;
format long;


%% 南方电网：广东电网分区11，贵州广西电网分区22，云南电网分区33
[system, bus, branch, Dbus, OLTC, Dbranch, gen, avr, pss, gov, HVDC, fault, Y]=readcase('case39');

%% 华中电网分区号为01，华东电网分区号为02，四川电网分区号为03
v = bus(:,2);
o = bus(:,3);
tic;
[bus] = Initial_Value(system, bus, branch, Dbus, Dbranch, Y, 1, 5);     %1-PQ，2-Gauss；3-Iteration Time
[bus, HVDC] = powerflow(system, bus, branch, Dbus, OLTC, Dbranch, HVDC, Y, gen, 2);      %1-PQ，2-NL
[savebus,savebranch,group]=readequ('DyEqucase');
[EQUbus,EQUbranch,EQUgen,EQUY,EQUfault] = EQU(bus,branch,Y,gen,savebus,savebranch,group,fault);

%% 等值后系统潮流计算
%[system, bus, branch, Dbus, OLTC, Dbranch, gen, avr, pss, gov, HVDC, fault, Y]=readcase('equdata');
%v=ones(size(bus,1),1);
%o=zeros(size(bus,1),1);
%[bus] = Initial_Value(system, bus, branch, Dbus, Dbranch, Y, 1, 5);     %1-PQ，2-Gauss；3-Iteration Time
%[bus, HVDC] = powerflow(system, bus, branch, Dbus, OLTC, Dbranch, HVDC, Y, gen, 2);       %1-PQ，2-NL
toc;

%[gen, HVDC, g, result] = tsp11(bus, gen, avr, pss, gov, fault, HVDC, Y, system);

%[PKE, PPE, PEF] = seekPKE(result, g, gen, system);

%plotdelta(result, gen, HVDC, fault, 0);