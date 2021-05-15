function [bus, HVDC, iter] = NL_powerflow(system, bus, Dbus, OLTC, Dbranch, HVDC, Y)
%% 电流接口的潮流程序（定功率运行，整流侧定电压，整流侧定触发角度，逆变侧定熄弧角）
% 注意，直流节点两侧的交流节点如果是PV节点的话，则无功部分没有注入元素，因此，最终求解结果交流侧是所得解，而直流侧各量并不一定是最终解；
% 即偏差判断时，F中包括了box中的P偏差，而没有包括Q偏差，因此有可能导致Q偏差不在允许偏差外。
threshold = 1e-10;                                      flags = 0;
j = sqrt(-1);                                           bus_num = size(bus, 1);
v = bus(:, 2);                                          o = bus(:, 3);
U = v.*(cos(o) + j*sin(o));                             S = bus(:, 4) - bus(:, 6) + j*(bus(:, 5) - bus(:, 7));
%% 直流系统初始化
if ~isempty(Dbus)
    %% 数据初始化
    Dln      = size(Dbranch, 1);                        unit     = ones(Dln, 1);
    rbus     = Dbranch(:,1);                            ibus     = Dbranch(:,2);
    racbus   = Dbus(rbus,2);                            iacbus   = Dbus(ibus,2);
    Nbr      = Dbus(rbus,4);                            Nbi      = Dbus(ibus,4);
    Vdropr   = Dbus(rbus,8)./(1000*HVDC.UdcB);          Vdropi   = Dbus(ibus,8)./(1000*HVDC.UdcB);
    Idn      = Dbranch(:,3)./(1000*HVDC.IdcB);          R        = Dbranch(:,4)./HVDC.ZdcB;
    L        = Dbranch(:,5)./(1000*HVDC.ZdcB);          C        = Dbranch(:,6).*HVDC.ZdcB;
    Ppos     = Dbranch(:,7);                            Pds      = Dbranch(:,8)./HVDC.MVABASE;
    Uds      = Dbranch(:,9)./HVDC.UdcB;
    alphar   = Dbranch(:,10)*pi/180;                    alphai   = ((180-20).*unit-Dbranch(:,11))*pi/180;
    gammar   = ((180-20).*unit-Dbranch(:,10))*pi/180;   gammai   = Dbranch(:,11)*pi/180;
    miur     = 20*pi/180.*unit;                         miui     = 20*pi/180.*unit;
    index    = find(Ppos == 2);
    if ~isempty(index)
        Idx        = (Uds(index)-sqrt(Uds(index).*Uds(index)-4.*R(index).*Pds(index)))./(2.*R(index));
        Pds(index) = Uds(index).*Idx;
    end
    %% 换流变压器折算，输入结构（--K:1--r+jx--）
    KarB   = Dbus(rbus,3)./HVDC.UdcB;                   KaiB   = Dbus(ibus,3)./HVDC.UdcB;
    OLTC(rbus,6) = OLTC(rbus,6)./KarB;                  OLTC(ibus,6) = OLTC(ibus,6)./KaiB;
    OLTC(rbus,7) = OLTC(rbus,7)./KarB;                  OLTC(ibus,7) = OLTC(ibus,7)./KaiB;
    OLTC(rbus,8) = OLTC(rbus,8)./KarB;                  OLTC(ibus,8) = OLTC(ibus,8)./KaiB;
    index = find(OLTC(:, 9) == 1);
    if ~isempty(index)
        OLTC(index, 6) = 1./OLTC(index, 6);             OLTC(index, 7) = 1./OLTC(index, 7);
        OLTC(index, 8) = 1./OLTC(index, 8);
    end
    index = find(OLTC(:, 9) == 2);
    if ~isempty(index)
        OLTC(index, 2) = OLTC(index, 2).*OLTC(index, 6).*OLTC(index, 6);
        OLTC(index, 3) = OLTC(index, 3).*OLTC(index, 6).*OLTC(index, 6);
    end
    Kars     = OLTC(rbus,6);                            Kais     = OLTC(ibus,6);
    Rcr      = OLTC(rbus,2).*KarB.^2;                   Rci      = OLTC(ibus,2).*KaiB.^2;
    Xcr      = OLTC(rbus,3).*KarB.^2;                   Xci      = OLTC(ibus,3).*KaiB.^2;

%% 直流系统初值设定
    Id       = Pds./Uds;                                Ids      = Id;
    Udr      = Uds;                                     Udi      = Udr - R.*Id;
    tanfir   = 0.5*unit;                                tanfii   = 0.5*unit;
    alpha_s  = alphar;                                  gamma_s  = gammai;
    Kar      = Kars;                                    Kai      = Kais;
    Udro     = 3*sqrt(2)/pi.*Kar.*Nbr.*bus(racbus,2);   Udio     = 3*sqrt(2)/pi.*Kai.*Nbi.*bus(iacbus,2);
    %% 矩阵位置初始化
    sn       = 19;
    Dindex   = zeros(Dln, 1);                           Dindex(1:Dln, 1)  = (1: sn: (sn*Dln))-1;
    Aindex   = zeros(Dln, 1);                           Aindex(1:Dln, 1)  = (1: 4: (4*Dln)) - 1;
    Irx   = zeros(Dln, 1);                              Iry      = zeros(Dln, 1);
    Iix   = zeros(Dln, 1);                              Iiy      = zeros(Dln, 1);
    DJI   = [
             Dindex+ 1*unit;    Dindex+ 1*unit;    Dindex+ 1*unit;    Dindex+ 1*unit;
             Dindex+ 2*unit;    Dindex+ 2*unit;    Dindex+ 2*unit;    Dindex+ 2*unit;    Dindex+ 2*unit;
             Dindex+ 3*unit;    Dindex+ 3*unit;    Dindex+ 3*unit;    Dindex+ 3*unit;
             Dindex+ 4*unit;    Dindex+ 4*unit;    Dindex+ 4*unit;    Dindex+ 4*unit;    Dindex+ 4*unit;
             Dindex+ 5*unit;    Dindex+ 5*unit;    Dindex+ 5*unit;
             Dindex+ 6*unit;    Dindex+ 6*unit;
             Dindex+ 7*unit;    Dindex+ 7*unit;    Dindex+ 7*unit;    Dindex+ 7*unit;
             Dindex+ 8*unit;    Dindex+ 8*unit;    Dindex+ 8*unit;    Dindex+ 8*unit;
             Dindex+ 9*unit;    Dindex+ 9*unit;    Dindex+ 9*unit;
             Dindex+10*unit;
             Dindex+11*unit;    Dindex+11*unit;    Dindex+11*unit;
             Dindex+12*unit;    Dindex+12*unit;
             Dindex+13*unit;    Dindex+13*unit;    Dindex+13*unit;    Dindex+13*unit;
             Dindex+14*unit;    Dindex+14*unit;    Dindex+14*unit;    Dindex+14*unit;
             Dindex+15*unit;    Dindex+15*unit;    Dindex+15*unit;
             Dindex+16*unit;    Dindex+16*unit;    Dindex+16*unit;
             Dindex+17*unit;
             Dindex+18*unit;
             Dindex+19*unit;    Dindex+19*unit;
        ];
    DJJ   = [
             Dindex+ 1*unit;    Dindex+ 2*unit;    Dindex+ 5*unit;    Dindex+ 8*unit;
             Dindex+ 1*unit;    Dindex+ 2*unit;    Dindex+ 5*unit;    Dindex+ 8*unit;    Dindex+10*unit;
             Dindex+ 3*unit;    Dindex+ 4*unit;    Dindex+ 5*unit;    Dindex+15*unit;
             Dindex+ 3*unit;    Dindex+ 4*unit;    Dindex+ 5*unit;    Dindex+15*unit;    Dindex+17*unit;
             Dindex+ 5*unit;    Dindex+ 8*unit;    Dindex+15*unit;
             Dindex+ 6*unit;    Dindex+ 7*unit;
             Dindex+ 5*unit;    Dindex+ 6*unit;    Dindex+ 8*unit;    Dindex+11*unit;
             Dindex+ 5*unit;    Dindex+ 6*unit;    Dindex+ 9*unit;    Dindex+11*unit;
             Dindex+ 9*unit;    Dindex+10*unit;    Dindex+11*unit;
             Dindex+11*unit;
             Dindex+ 9*unit;    Dindex+11*unit;    Dindex+12*unit;
             Dindex+13*unit;    Dindex+14*unit;
             Dindex+ 5*unit;    Dindex+13*unit;    Dindex+15*unit;    Dindex+18*unit;
             Dindex+ 5*unit;    Dindex+13*unit;    Dindex+16*unit;    Dindex+18*unit;
             Dindex+16*unit;    Dindex+17*unit;    Dindex+19*unit;
             Dindex+16*unit;    Dindex+18*unit;    Dindex+19*unit;
             Dindex+19*unit;
             Dindex+ 8*unit;
             Dindex+ 5*unit;    Dindex+ 8*unit;
        ];
     DC2ACI= [
              Dindex+ 1*unit;          Dindex+ 1*unit;
              Dindex+ 2*unit;          Dindex+ 2*unit;
              Dindex+ 3*unit;          Dindex+ 3*unit;
              Dindex+ 4*unit;          Dindex+ 4*unit;
              Dindex+ 6*unit;          Dindex+ 6*unit;
              Dindex+12*unit;          Dindex+12*unit;
        ];
    DC2ACJ= [
              Aindex+ 1*unit;         Aindex+ 2*unit;
              Aindex+ 1*unit;         Aindex+ 2*unit;
              Aindex+ 3*unit;         Aindex+ 4*unit;
              Aindex+ 3*unit;         Aindex+ 4*unit;
              Aindex+ 1*unit;         Aindex+ 2*unit;
              Aindex+ 3*unit;         Aindex+ 4*unit;
        ];
    DFI   = [
             Dindex+ 1*unit;
             Dindex+ 2*unit;
             Dindex+ 3*unit;
             Dindex+ 4*unit;
             Dindex+ 5*unit;
             Dindex+ 6*unit;
             Dindex+ 7*unit;
             Dindex+ 8*unit;
             Dindex+ 9*unit;
             Dindex+10*unit;
             Dindex+11*unit;
             Dindex+12*unit;
             Dindex+13*unit;
             Dindex+14*unit;
             Dindex+15*unit;
             Dindex+16*unit;
             Dindex+17*unit;
             Dindex+18*unit;
             Dindex+19*unit;
        ];
    DFJ   = ones(sn*Dln, 1);
end
%% 潮流计算
while flags == 0
    pqbus = find(bus(:,13) == 1 );                      pqn = size(pqbus,1);
    pvbus = find(bus(:,13) == 2 );                      pvn = size(pvbus,1);
    i1 = 1:pqn;                                         i2 = (pqn+1):2*pqn;
    i3 = (2*pqn+1):(2*pqn+pvn);                         i4 = (2*pqn+pvn+1):(2*pqn+2*pvn);
    F  = 1;                                             iter = 0;
    while max(abs(F)) > threshold
       %% 交流系统
        Ux   = real(U);                                 Uy   = imag(U);
        diagUx = sparse(1:bus_num,1:bus_num,Ux,bus_num,bus_num);
        diagUy = sparse(1:bus_num,1:bus_num,Uy,bus_num,bus_num);
        ReY = real(Y);                                  ImY  = imag(Y);
        YU  = Y*U;                                      zmatrix = sparse([],[],[],pvn,pqn);
        P2Ux = -diagUx*ReY-diagUy*ImY+sparse(1:bus_num,1:bus_num,-real(YU),bus_num,bus_num);
        P2Uy =  diagUx*ImY-diagUy*ReY+sparse(1:bus_num,1:bus_num,-imag(YU),bus_num,bus_num);
        Q2Ux =  diagUx*ImY-diagUy*ReY+sparse(1:bus_num,1:bus_num, imag(YU),bus_num,bus_num);
        Q2Uy =  diagUx*ReY+diagUy*ImY+sparse(1:bus_num,1:bus_num,-real(YU),bus_num,bus_num);
        U22Ux= sparse(1:bus_num,1:bus_num,-2*Ux,bus_num,bus_num);
        U22Uy= sparse(1:bus_num,1:bus_num,-2*Uy,bus_num,bus_num);
        misS = S-U.*conj(Y*U);                          misV = v.*v-(Ux.*Ux+Uy.*Uy);
        if ~isempty(Dbus)
           %% 直流系统
            Ix2Ux    = sparse([],[],[],bus_num,bus_num);    Ix2Uy    = sparse([],[],[],bus_num,bus_num);
            Iy2Ux    = sparse([],[],[],bus_num,bus_num);    Iy2Uy    = sparse([],[],[],bus_num,bus_num);
            Ix       = zeros(bus_num, 1);                   Iy       = zeros(bus_num, 1);
            CIx      = zeros(bus_num, 1);                   CIy      = zeros(bus_num, 1);
            Urx      = Ux(racbus);                          Ury      = Uy(racbus);
            Uix      = Ux(iacbus);                          Uiy      = Uy(iacbus);
            Ur       = sqrt(Urx.*Urx+Ury.*Ury);             Ui       = sqrt(Uix.*Uix+Uiy.*Uiy);
            
            % 更新直流系统对交流侧电流的注入
            for ii = 1:Dln
                i = racbus(ii);
                Ix(i) = Ix(i) + Irx(ii);                    Iy(i) = Iy(i) + Iry(ii);
            end
            for ii = 1:Dln
                i = iacbus(ii);
                Ix(i) = Ix(i) + Iix(ii);                    Iy(i) = Iy(i) + Iiy(ii);
            end
            diagIx = sparse(1:bus_num,1:bus_num,Ix,bus_num,bus_num);
            diagIy = sparse(1:bus_num,1:bus_num,Iy,bus_num,bus_num);
            I      = Ix + j.*Iy;
            DJV    = [
                      Urx;    Ury;    Udr;    Id;
                      Ury;    -Urx;   Udr.*tanfir;    Id.*tanfir;    Udr.*Id;
                      Uix;    Uiy;    -Udi;   -Id;
                      Uiy;    -Uix;   Udi.*tanfii;    Id.*tanfii;    Udi.*Id;
                      -R;     unit;   -unit;
                      unit;   -3*sqrt(2)/pi.*Nbr.*Ur;
                      Nbr.^2.*(3/pi.*Xcr+2*Rcr);    -cos(alphar);    unit;    Udro.*sin(alphar);
                      6/pi.*Nbr.^2.*Xcr;    cos(alphar+miur)-cos(alphar);    -Udro.*sin(alphar+miur);    Udro.*(-sin(alphar+miur)+sin(alphar));
                      2.*tanfir.*sin(2*alphar+2*miur)-2.*unit+2*cos(2*alphar+2*miur);    cos(2*alphar)-cos(2*alphar+2*miur);    tanfir.*(-2*sin(2*alphar)+2*sin(2*alphar+2*miur))-2*cos(2*alphar)+2*cos(2*alphar+2*miur);
                      unit;
                      unit;    unit;    unit;
                      unit;    -3*sqrt(2)/pi.*Nbi.*Ui;
                      -Nbi.^2.*(3/pi.*Xci+2*Rci);    cos(alphai);    unit;    -Udio.*sin(alphai);
                      6/pi.*Nbi.^2.*Xci;    cos(alphai+miui)-cos(alphai);    -Udio.*sin(alphai+miui);    Udio.*(-sin(alphai+miui)+sin(alphai));
                      2.*tanfii.*sin(2*gammai+2*miui)-2.*unit+2*cos(2*gammai+2*miui);    cos(2*gammai)-cos(2*gammai+2*miui);    tanfii.*(-2*sin(2*gammai)+2*sin(2*gammai+2*miui))-2*cos(2*gammai)+2*cos(2*gammai+2*miui);
                      unit;    unit;    unit;
                      unit;
                      unit;
                      -Udr;    -Id;
                ];
            DC2ACV = [
                      Irx;        Iry;
                      -Iry;       Irx;
                      Iix;        Iiy;
                      -Iiy;       Iix;
                      -3*sqrt(2)/pi.*Kar.*Nbr.*Urx./Ur;    -3*sqrt(2)/pi.*Kar.*Nbr.*Ury./Ur;
                      -3*sqrt(2)/pi.*Kai.*Nbi.*Uix./Ui;    -3*sqrt(2)/pi.*Kai.*Nbi.*Uiy./Ui;
                ];
            DFV    = -[
                        Urx.*Irx + Ury.*Iry + Udr.*Id;
                        Ury.*Irx - Urx.*Iry + Udr.*Id.*tanfir;
                        Uix.*Iix + Uiy.*Iiy - Udi.*Id;
                        Uiy.*Iix - Uix.*Iiy + Udi.*Id.*tanfii;
                        Udr - Udi - R.*Id;
                        Udro - 3*sqrt(2)/pi.*Kar.*Nbr.*Ur;
                        Udr - Udro.*cos(alphar) + Nbr.^2.*(3/pi*Xcr+2*Rcr).*Id + 2*Nbr.*Vdropr;
                        Udro.*cos(alphar+miur) - Udro.*cos(alphar) + 6/pi.*Xcr.*Nbr.^2.*Id;
                        tanfir.*(cos(2*alphar)-cos(2*alphar+2*miur)) - (2*miur+sin(2*alphar)-sin(2*alphar+2*miur));
                        alphar - alpha_s;
                        alphar + miur + gammar - pi;
                        Udio - 3*sqrt(2)/pi.*Kai.*Nbi.*Ui;
                        Udi + Udio.*cos(alphai) - Nbi.^2.*(3/pi*Xci+2*Rci).*Id - 2*Nbi.*Vdropi;
                        Udio.*cos(alphai+miui) - Udio.*cos(alphai) + 6/pi.*Xci.*Nbi.^2.*Id;
                        tanfii.*(cos(2*gammai)-cos(2*gammai+2*miui)) - (2*miui+sin(2*gammai)-sin(2*gammai+2*miui));
                        alphai + miui + gammai - pi;
                        gammai - gamma_s;
                        Udr - Uds;
                        Pds - Udr.*Id;
                ];
            DJ    = sparse(DJI, DJJ, DJV, sn*Dln, sn*Dln);
            DC2AC = sparse(DC2ACI, DC2ACJ, DC2ACV, sn*Dln, 4*Dln);
            DF    = sparse(DFI, DFJ, DFV, sn*Dln, 1);
            DK    = -DJ\DC2AC;
            DC    = DJ\DF;
            %% 解决换流节点上所连直流母线的个数大于1的情况
            for ii = 1:Dln
                i  = racbus(ii);
                for jj = 1:Dln
                    j = racbus(jj);
                    Ix2Ux(i, j) = Ix2Ux(i, j) + DK(Dindex(ii)+1, Aindex(jj)+1);
                    Ix2Uy(i, j) = Ix2Uy(i, j) + DK(Dindex(ii)+1, Aindex(jj)+2);
                    Iy2Ux(i, j) = Iy2Ux(i, j) + DK(Dindex(ii)+2, Aindex(jj)+1);
                    Iy2Uy(i, j) = Iy2Uy(i, j) + DK(Dindex(ii)+2, Aindex(jj)+2);
                end
            end
            for ii = 1:Dln
                i  = iacbus(ii);
                for jj = 1:Dln
                    j = iacbus(jj);
                    Ix2Ux(i, j) = Ix2Ux(i, j) + DK(Dindex(ii)+3, Aindex(jj)+3);
                    Ix2Uy(i, j) = Ix2Uy(i, j) + DK(Dindex(ii)+3, Aindex(jj)+4);
                    Iy2Ux(i, j) = Iy2Ux(i, j) + DK(Dindex(ii)+4, Aindex(jj)+3);
                    Iy2Uy(i, j) = Iy2Uy(i, j) + DK(Dindex(ii)+4, Aindex(jj)+4);
                end
            end
            for ii = 1:Dln
                i  = racbus(ii);
                for jj = 1:Dln
                    j = iacbus(jj);
                    Ix2Ux(i, j) = Ix2Ux(i, j) + DK(Dindex(ii)+1, Aindex(jj)+3);
                    Ix2Uy(i, j) = Ix2Uy(i, j) + DK(Dindex(ii)+1, Aindex(jj)+4);
                    Iy2Ux(i, j) = Iy2Ux(i, j) + DK(Dindex(ii)+2, Aindex(jj)+3);
                    Iy2Uy(i, j) = Iy2Uy(i, j) + DK(Dindex(ii)+2, Aindex(jj)+4);
                end
            end
            for ii = 1:Dln
                i  = iacbus(ii);
                for jj = 1:Dln
                    j = racbus(jj);
                    Ix2Ux(i, j) = Ix2Ux(i, j) + DK(Dindex(ii)+3, Aindex(jj)+1);
                    Ix2Uy(i, j) = Ix2Uy(i, j) + DK(Dindex(ii)+3, Aindex(jj)+2);
                    Iy2Ux(i, j) = Iy2Ux(i, j) + DK(Dindex(ii)+4, Aindex(jj)+1);
                    Iy2Uy(i, j) = Iy2Uy(i, j) + DK(Dindex(ii)+4, Aindex(jj)+2);
                end
            end
            for ii = 1:Dln
                i = racbus(ii);
                CIx(i) = CIx(i) + DC(Dindex(ii)+1);            CIy(i) = CIy(i) + DC(Dindex(ii)+2);
            end
            for ii = 1:Dln
                i = iacbus(ii);
                CIx(i) = CIx(i) + DC(Dindex(ii)+3);            CIy(i) = CIy(i) + DC(Dindex(ii)+4);
            end
            %修改雅可比矩阵和偏差量元素
            KP2Ux = diagUx*Ix2Ux + diagUy*Iy2Ux;                KP2Uy = diagUx*Ix2Uy + diagUy*Iy2Uy;
            KQ2Ux = diagUy*Ix2Ux - diagUx*Iy2Ux;                KQ2Uy = diagUy*Ix2Uy - diagUx*Iy2Uy;
            CP    = diagUx*CIx + diagUy*CIy;                    CQ    = diagUy*CIx - diagUx*CIy;
            
            P2Ux = P2Ux + KP2Ux + diagIx;
            P2Uy = P2Uy + KP2Uy + diagIy;
            Q2Ux = Q2Ux + KQ2Ux - diagIy;
            Q2Uy = Q2Uy + KQ2Uy + diagIx;
            
            j    = sqrt(-1);
            misS = misS + U.*conj(I) + (CP +j*CQ);
        end
        %% 雅可比矩阵
        J    = [
                P2Ux(pqbus,pqbus)   P2Uy(pqbus,pqbus)   P2Ux(pqbus,pvbus)   P2Uy(pqbus,pvbus);
                Q2Ux(pqbus,pqbus)   Q2Uy(pqbus,pqbus)   Q2Ux(pqbus,pvbus)   Q2Uy(pqbus,pvbus);
                P2Ux(pvbus,pqbus)   P2Uy(pvbus,pqbus)   P2Ux(pvbus,pvbus)   P2Uy(pvbus,pvbus);
                zmatrix             zmatrix             U22Ux(pvbus,pvbus)  U22Uy(pvbus,pvbus)
                ];
        F    = -[
                real(misS(pqbus,1));
                imag(misS(pqbus,1));
                real(misS(pvbus,1));
                misV(pvbus,1)
                ];
        x   = J\F;
        if ~isempty(Dbus)
            dUxy = zeros(4*Dln, 1);
            dUx  = zeros(bus_num,1);                    dUy  = zeros(bus_num,1);
            dUx(pqbus) = x(i1);                         dUy(pqbus) = x(i2);
            dUx(pvbus) = x(i3);                         dUy(pvbus) = x(i4);
            dUrx = dUx(racbus);                         dUry = dUy(racbus);
            dUix = dUx(iacbus);                         dUiy = dUy(iacbus);
            dUxy(Aindex+1*unit, 1) = dUrx;              dUxy(Aindex+2*unit, 1) = dUry;
            dUxy(Aindex+3*unit, 1) = dUix;              dUxy(Aindex+4*unit, 1) = dUiy;
            dDx  = DK*dUxy + DC;
            
            Irx     = Irx    + dDx(Dindex+ 1*unit);     Iry     = Iry    + dDx(Dindex+ 2*unit);
            Iix     = Iix    + dDx(Dindex+ 3*unit);     Iiy     = Iiy    + dDx(Dindex+ 4*unit);
            Id      = Id     + dDx(Dindex+ 5*unit);
            Udro    = Udro   + dDx(Dindex+ 6*unit);     Udio    = Udio   + dDx(Dindex+13*unit);
            Kar     = Kar    + dDx(Dindex+ 7*unit);     Kai     = Kai    + dDx(Dindex+14*unit);
            Udr     = Udr    + dDx(Dindex+ 8*unit);     Udi     = Udi    + dDx(Dindex+15*unit);
            miur    = miur   + dDx(Dindex+ 9*unit);     miui    = miui   + dDx(Dindex+16*unit);
            tanfir  = tanfir + dDx(Dindex+10*unit);     tanfii  = tanfii + dDx(Dindex+17*unit);
            alphar  = alphar + dDx(Dindex+11*unit);     alphai  = alphai + dDx(Dindex+18*unit);
            gammar  = gammar + dDx(Dindex+12*unit);     gammai  = gammai + dDx(Dindex+19*unit);
        end
        %修正变量
        Ux(pqbus) = Ux(pqbus) + x(i1);                  Uy(pqbus) = Uy(pqbus) + x(i2);
        Ux(pvbus) = Ux(pvbus) + x(i3);                  Uy(pvbus) = Uy(pvbus) + x(i4);
        U    = Ux + j*Uy;
        iter = iter + 1;
        if iter > 300
            error('The Newton-Raphson power flow does not converge!');
        end
    end
    if ~isempty(Dbus)
        Kar_min = OLTC(rbus, 7);                        Kar_max = OLTC(rbus, 8);
        Kai_min = OLTC(ibus, 7);                        Kai_max = OLTC(ibus, 8);
        index1  = find(Kar>Kar_max);                    index2  = find(Kai>Kai_max);
        index3  = find(Kar<Kar_min);                    index4  = find(Kai<Kai_min);
        if ~isempty(index1)
            Kars(index1,1) = Kar_max(index1,1);
            fprintf('The rectifier transformer ratio is upper than the maximun limit!\n');
            fprintf('%5d \n', index1');
        end
        if ~isempty(index2)
            Kais(index2,1) = Kai_max(index2,1);
            fprintf('The inverter transformer ratio is upper than the maximun limit!\n');
            fprintf('%5d \n', index2');
        end
        if ~isempty(index3)
            Kars(index3,1) = Kar_min(index3,1);
            fprintf('The rectifier transformer ratio is lower than the minimun limit!\n');
            fprintf('%5d \n', index3');
        end
        if ~isempty(index4)
            Kais(index4,1) = Kai_min(index4,1);
            fprintf('The inverter transformer ratio is lower than the minimun limit!\n');
            fprintf('%5d \n', index4');
        end
    end
    flags = 1;
    Qe=imag(U.*conj(Y*U));                         gq   =Qe+bus(:,7);
    gqmax=bus(:,8);                                gqmin=bus(:,9);
    index1=find(gq>gqmax);                         index2=find(bus(:,13)==2);          index=intersect(index1,index2);
    if size(index,1)==0
        flags=1;
    else
        flags=0;
        bus(index,13)=1;
        bus(index,5) = bus(index,8);
        S = (bus(:,4) - bus(:,6))+j*(bus(:,5)- bus(:,7));
        pqbus	= find(bus(:,13) == 1 );              pvbus	= find(bus(:,13) == 2 );
    end
    index1=find(gq<gqmin);                         index=intersect(index1,index2);
    if size(index,1)==0
        flags=1;
    else
        flags=0;
        bus(index,13)=1;
        bus(index,5) = bus(index,9);
        S = (bus(:,4) - bus(:,6))+j*(bus(:,5)- bus(:,7));
        pqbus	= find(bus(:,13) == 1 );              pvbus	= find(bus(:,13) == 2 );
    end
    fprintf('The Newton-Raphson power flow converge, iter = %2d\n',iter);
end
%% 输出显示  17-发电机总的输出功率
v  = abs(U);                                            o = angle(U);
bus(:,2) = v;                                           bus(:,3) = o;
Se = U.*conj(Y*U);
bus(:,16) = real(Se);                                   bus(:,17) = imag(Se);
bus(:,18) = Se + bus(:,6) + j*bus(:,7);

%% 直流系统计算结果，标幺值
if ~isempty(Dbus)
    HVDC.rbus   = rbus;                                 HVDC.ibus   = ibus;
    HVDC.racbus = racbus;                               HVDC.iacbus = iacbus;
    HVDC.Nbr    = Nbr;                                  HVDC.Nbi    = Nbi;
    HVDC.Rcr    = Rcr;                                  HVDC.Rci    = Rci;
    HVDC.Xcr    = Xcr;                                  HVDC.Xci    = Xci;
    HVDC.Vdropr = Vdropr;                               HVDC.Vdropi = Vdropi;
    HVDC.Pds    = Pds;                                  HVDC.Uds    = Uds;
    
    HVDC.Idr    = Id;                                   HVDC.Idi    = Id;
    HVDC.Udr    = Udr;                                  HVDC.Udi    = Udi;
    HVDC.Udro   = Udro;                                 HVDC.Udio   = Udio;
    HVDC.Kar    = Kar;                                  HVDC.Kai    = Kai;
    HVDC.tanfir = tanfir;                               HVDC.tanfii = tanfii;
    HVDC.alphar = alphar;                               HVDC.alphai = alphai;
    HVDC.miur   = miur;                                 HVDC.miui   = miui;
    HVDC.gammar = gammar;                               HVDC.gammai = gammai;
    HVDC.R      = R;                                    HVDC.L      = L;
    HVDC.C      = C;                                    HVDC.Idn    = Idn;
    
    HVDC.Irx    = Irx;                                  HVDC.Iry    = Iry;
    HVDC.Iix    = Iix;                                  HVDC.Iiy    = Iiy;
    
    HVDC.Lsr    = Dbus(rbus,5)./1000;                   HVDC.Lsi    = Dbus(ibus,5)./1000;
    HVDC.Aminr  = Dbus(rbus,6)*pi/180;                  HVDC.Amini  = Dbus(ibus,6)*pi/180;
    HVDC.Astopr = Dbus(rbus,7)*pi/180;                  HVDC.Astopi = Dbus(ibus,7)*pi/180;
    HVDC.Uc     = Udr-0.5*R.*Id;                        HVDC.Ids    = Ids;
    HVDC.Gimin  = 16*pi/180*unit;                       HVDC.Gimax  = 60*pi/180*unit;
    HVDC.Armin  = 5*pi/180*unit;                        HVDC.Armax  = 90*pi/180*unit;
    HVDC.krr    = 0.5*(cos(alphar)+cos(alphar+miur)).*sqrt(1+(miur.*csc(miur).*csc(2*alphar+miur)-cot(2*alphar+miur)).^2);
    HVDC.kri    = 0.5*(cos(gammai)+cos(gammai+miui)).*sqrt(1+(miui.*csc(miui).*csc(2*gammai+miui)-cot(2*gammai+miui)).^2);
    HVDC.L      = HVDC.L + HVDC.Lsr./HVDC.ZdcB + HVDC.Lsi./HVDC.ZdcB;
    HVDC.C      = zeros(Dln,1);                         HVDC.Dln    = Dln;
    HVDC.IsBlock= zeros(Dln,1);
end
return;