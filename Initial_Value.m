function [bus] = Initial_Value(system, bus, branch, Dbus, Dbranch, Y, method, iter)
%% PQ�ֽⷨ���߸�˹���¶��������趨��������������ֵ����������ƽ�����������������⣩
Pgen    = bus(:, 4);                                            Qgen    = bus(:, 5);
Pload   = bus(:, 6);                                            Qload   = bus(:, 7);
%% %ֱ��ϵͳ�г�����Ч�ɸ��ɺͷ���
if ~isempty(Dbus)
    MVABASE = system.MVABASE;
    bus_id   = Dbus(:, 2);                                      Dln     = size(Dbranch, 1);
    ino      = Dbranch(:, 1);                                   jno     = Dbranch(:, 2);
    R_LD     = Dbranch(:, 4);                                   D_set   = Dbranch(:, 7);
    Pds      = Dbranch(:, 8);                                   Uds     = Dbranch(:, 9);
    
  	V_RAB    = Dbus(ino, 10);                                   V_IAB   = Dbus(jno, 10);
    VB1      = V_RAB;                                           VB2     = V_IAB;
    Kv       = VB2./VB1;
    IB       = MVABASE./VB1;                                    ZB      = VB1./IB;
    
    Pds      = Pds/MVABASE;                                     Uds     = Uds./VB1;
    R_LD     = R_LD./ZB;                                        Pdr     = Pds;
    
    index    = find(D_set == 2);    %������������
    if ~isempty(index)
        Vdi  = (Uds + sqrt(Uds.*Uds - 4.*Pds.*R_LD))./(2.*Kv);  Id      = (Uds - Kv.*Vdi)./R_LD;
        Ptpr = Uds.*Id;                                         Pdr(index) = Ptpr(index);
    end
    Vdi     = (Uds.*Uds - Pdr.*R_LD)./(Kv.*Uds);                Id      = (Uds - Kv.*Vdi)./R_LD;
  	Pdi     = Kv.*Vdi.*Id;
    Qdr     = 0.5*Pdr;                                       	Qdi     = 0.5*Pdi;
    
    for i = 1:Dln
        fno = bus_id(ino(i));                                   tno     = bus_id(jno(i));
        Pload(fno) = Pload(fno) + Pdr(i);                       Qload(fno) = Qload(fno) + Qdr(i);
        Pgen(tno)  = Pgen(tno)  + Pdi(i);                       Qgen(tno)  = Qgen(tno)  + Qdi(i);
    end
end
j       = sqrt(-1);
S       = (Pgen - Pload) + j*(Qgen - Qload);
%% ��ֵ�����趨
if method == 1
    [bus] = Initial_PQ_powerflow(S, bus, branch, Y, iter);
elseif method == 2
    [bus] = Initial_GaussSiedel(S, bus, Y, iter);
else
    error('The method is not true!!!');
end
%% PQ�ֽⷨ
function [bus] = Initial_PQ_powerflow(S, bus, branch, Y, iter)
j       = sqrt(-1);
inP     = real(S);                                      inQ     = imag(S);
Um      = bus(:, 2);                                   	Ua      = bus(:, 3);
refbus  = find(bus(:, 13) == 3);                        pvbus   = find(bus(:, 13) == 2);
nl      = size(branch, 1);                            	nb      = size(bus, 1);
f       = branch(:, 1);                                	t       = branch(:, 2);
Cf      = sparse(f, 1:nl, ones(nl, 1), nb, nl);       	Ct      = sparse(t, 1:nl, ones(nl, 1), nb, nl);
gl      = find(f == t);
Cf(:, gl) = 0;                                          Ct(:, gl) = 0;
X       = branch(:, 4);                               	X(X == 0) = 1e-10;
Bs      = 1./X;
Btt     = Bs;                                       	Bff     = Bs;
Bft     = -Bs;                                      	Btf     = -Bs;
B       = Cf * spdiags(Bff, 0, nl, nl) * Cf' + Cf * spdiags(Bft, 0, nl, nl) * Ct' + ...
          Ct * spdiags(Btf, 0, nl, nl) * Cf' + Ct * spdiags(Btt, 0, nl, nl) * Ct';
BB      = -imag(Y);
%�ſɱȾ���
Jp      = B;                                            Jq      = BB;
Jp(refbus,:) = 0;                                       Jq(refbus,:) = 0;
Jq(pvbus,:)  = 0;
dJvv    = sparse(refbus, refbus, 1, nb, nb);        	dJpv    = sparse(pvbus, pvbus, 1, nb, nb);
Jp      = Jp + dJvv;                                    Jq      = Jq + dJvv + dJpv;
%������ʼ
K       = 0;
while 1
    %�й�
	U   = Um.*(cos(Ua) + j*sin(Ua));                  	Pe  = real(U.*conj(Y*U));
  	dP  = (inP - Pe)./Um;                             	dP(refbus) = 0;
   	dUa = Jp\dP;                                      	Ua  = Ua + dUa;
    %�޹�
   	U   = Um.*(cos(Ua) + j*sin(Ua));                  	Qe  = imag(U.*conj(Y*U));
  	dQ  = (inQ - Qe)./Um;                             	dQ(refbus) = 0;
	dQ(pvbus) = 0;
  	dUm = Jq\dQ;                                        Um  = Um + dUm;
    %������������
    K   = K + 1;
    if K > iter
        break;
    end
end
bus(:, 2) = Um;                                         bus(:, 3) = Ua;
%% ��˹���¶�������
function [bus] = Initial_GaussSiedel(S, bus, Y, iter)
% PQ�ֽⷨ�Ѿ��������������������������ֵ���㣬�ʸ�˹���¶���������δ��ӣ�����Ĵ���δ��֤��
% ���ɾ�����ʽ�ļ����ٶȽϸߣ��������Լ���迹������ʽ�ļ����ٶȼ��ͣ��������ԽϺ�
bus_num = size(bus, 1);                                 j       = sqrt(-1);
Um      = bus(:, 2);                                   	Ua      = bus(:, 3);
U       = Um.*(cos(Ua) + j*sin(Ua));                    Z       = inv(Y);
Se      = S;
%������ʼ
K       = 0;
while 1
    for i = 1:bus_num
        if bus(i, 13) == 3
            %ƽ��ڵ㣬������
        elseif bus(i, 13) == 2
            diagU = sparse(1:bus_num, 1:bus_num, 1./conj(U), bus_num, bus_num);
            U(i) = Z(i, :) * diagU * conj(Se) - Z(i, i) * conj(Se(i))/conj(U(i));
            U(i) = Um(i) * U(i)/abs(U(i));
        else
            diagU = sparse(1:bus_num, 1:bus_num, 1./conj(U), bus_num, bus_num);
            U(i) = Z(i, :) * diagU * conj(Se) - Z(i, i) * conj(Se(i))/conj(U(i));
        end
    end
    %������������
    K   = K + 1;
    if K > iter
        break;
    end
end
bus(:, 2) = abs(U);                                   	bus(:, 3) = angle(U);
%% ����
return;