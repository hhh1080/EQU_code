clear;
clc;
close all;
format long;


%% �Ϸ��������㶫��������11�����ݹ�����������22�����ϵ�������33
[system, bus, branch, Dbus, OLTC, Dbranch, gen, avr, pss, gov, HVDC, fault, Y]=readcase('case39');

%% ���е���������Ϊ01����������������Ϊ02���Ĵ�����������Ϊ03
v = bus(:,2);
o = bus(:,3);
tic;
[bus] = Initial_Value(system, bus, branch, Dbus, Dbranch, Y, 1, 5);     %1-PQ��2-Gauss��3-Iteration Time
[bus, HVDC] = powerflow(system, bus, branch, Dbus, OLTC, Dbranch, HVDC, Y, gen, 2);      %1-PQ��2-NL
[savebus,savebranch,group]=readequ('DyEqucase');
[EQUbus,EQUbranch,EQUgen,EQUY,EQUfault] = EQU(bus,branch,Y,gen,savebus,savebranch,group,fault);

%% ��ֵ��ϵͳ��������
%[system, bus, branch, Dbus, OLTC, Dbranch, gen, avr, pss, gov, HVDC, fault, Y]=readcase('equdata');
%v=ones(size(bus,1),1);
%o=zeros(size(bus,1),1);
%[bus] = Initial_Value(system, bus, branch, Dbus, Dbranch, Y, 1, 5);     %1-PQ��2-Gauss��3-Iteration Time
%[bus, HVDC] = powerflow(system, bus, branch, Dbus, OLTC, Dbranch, HVDC, Y, gen, 2);       %1-PQ��2-NL
toc;

%[gen, HVDC, g, result] = tsp11(bus, gen, avr, pss, gov, fault, HVDC, Y, system);

%[PKE, PPE, PEF] = seekPKE(result, g, gen, system);

%plotdelta(result, gen, HVDC, fault, 0);