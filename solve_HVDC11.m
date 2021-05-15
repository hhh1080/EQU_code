function [HVDC] = solve_HVDC11(system, faultflag, HVDC, IsLimited, times)
%% 交直流联立求解法
if ~isempty(HVDC)
    fhz     = system.fhz;              hh     = system.hh;
    dc_num = size(HVDC.ibus,1);
    HVDC.box   = zeros(HVDC.sp(dc_num,1)+HVDC.tsn(dc_num,1),1);
    HVDC.g4n   = zeros(HVDC.sp(dc_num,1)+HVDC.tsn(dc_num,1),4);
    dcyI = [];        dcyJ = [];        dcyV = [];
    
    %% 判断直流线路是否闭锁，VDCOL控制，整流侧电压低于0.9---定电流控制，整流侧电压高于0.95---定功率控制
    for i = 1:dc_num
        Ur   = sqrt(HVDC.Urx(i,1)*HVDC.Urx(i,1) + HVDC.Ury(i,1)*HVDC.Ury(i,1));
        Ui   = sqrt(HVDC.Uix(i,1)*HVDC.Uix(i,1) + HVDC.Uiy(i,1)*HVDC.Uiy(i,1));
        IsBlock = HVDC.IsBlock(i,1);      %闭锁标志位
        if ((Ur == 0)||(Ui == 0))
            HVDC.IsBlock(i,1) = 1;
            HVDC.newv(i,:) = 0;
            continue;
        else
            Udro = 3*sqrt(2)/pi*HVDC.Nbr(i,1).*HVDC.Kar(i,1).*Ur;
            Udio = 3*sqrt(2)/pi*HVDC.Nbi(i,1).*HVDC.Kai(i,1).*Ui;
            Udr  = Udro.*cos(HVDC.alphar(i,1)) - HVDC.Nbr(i,1).^2.*(3/pi*HVDC.Xcr(i,1)+2*HVDC.Rcr(i,1)).*HVDC.Idr(i,1) - 2*HVDC.Nbr(i,1).*HVDC.Vdropr(i,1);
            Udi  =-Udio.*cos(HVDC.alphai(i,1)) + HVDC.Nbi(i,1).^2.*(3/pi*HVDC.Xci(i,1)+2*HVDC.Rci(i,1)).*HVDC.Idi(i,1) + 2*HVDC.Nbi(i,1).*HVDC.Vdropi(i,1);
            minUd   = 0.1;
            if ((Udr < minUd)||(Udi < minUd))
                HVDC.IsBlock(i,1) = 1;
            else
                HVDC.IsBlock(i,1) = 0;
            end
        end
        
        if (IsBlock == 1)&&(HVDC.IsBlock(i,1) == 1)
            HVDC.newv(i,:) = 0;
        elseif (IsBlock == 1)&&(HVDC.IsBlock(i,1) == 0)
            HVDC.newv(i,:) = HVDC.Start(i,:);
        elseif (IsBlock == 0)&&(HVDC.IsBlock(i,1) == 1)
            HVDC.newv(i,:) = 0;
        elseif (IsBlock == 0)&&(HVDC.IsBlock(i,1) == 0)
        end
        
        Idr = HVDC.newv(i,5);
        Idi = HVDC.newv(i,6);
        if (Udr > 0.9 && Udr > Udi)
            Iord(i,1) = HVDC.Pds(i,1)./Udr;
            continue;
        end
        Ucomp = Udr - HVDC.R(i,1)/2*Idr;
         
        if (Ucomp > 0.9)
            Iord(i,1) = (Ucomp+0.1)*HVDC.Ids(i,1);
        elseif ((Ucomp <= 0.9)&&(Ucomp > 0.4))
            Iord(i,1) = (0.9*Ucomp+0.19)*HVDC.Ids(i,1);
        elseif ((Ucomp <= 0.4)&&(Ucomp > 0.2))
            Iord(i,1) = 0.55*HVDC.Ids(i,1);
        else
            Iord(i,1) = 0;
        end
    end
    
%     HVDC.IsBlock = [0;1;0];
    %% 直流仿真参数重新定义
    IsOn = find(HVDC.IsBlock == 0);
    if isempty(IsOn)
        fprintf('no HVDC works\n');
        return;
    end
    unit   = ones(size(IsOn,1),1);
    Nbr  = HVDC.Nbr(IsOn,1);     Nbi  = HVDC.Nbi(IsOn,1);       Iord   = Iord(IsOn,1);
    R    = HVDC.R(IsOn,1);       L    = HVDC.L(IsOn,1);         C      = HVDC.C(IsOn,1);
    Rcr  = HVDC.Rcr(IsOn,1);     Rci  = HVDC.Rci(IsOn,1);       Vdropr = HVDC.Vdropr(IsOn,1);
    Xcr  = HVDC.Xcr(IsOn,1);     Xci  = HVDC.Xci(IsOn,1);       Vdropi = HVDC.Vdropi(IsOn,1);
    Kars = HVDC.Kar(IsOn,1);     Kais = HVDC.Kai(IsOn,1);
    Kpr  = HVDC.Kpr(IsOn,1);     Kir  = HVDC.Kir(IsOn,1);
    Kpi  = HVDC.Kpi(IsOn,1);     Kii  = HVDC.Kii(IsOn,1);
    
    %% 
    newIrx    =   HVDC.newv(IsOn, 1);          oldIrx    =   HVDC.state(IsOn, 1);
    newIry    =   HVDC.newv(IsOn, 2);          oldIry    =   HVDC.state(IsOn, 2);
    newIix    =   HVDC.newv(IsOn, 3);          oldIix    =   HVDC.state(IsOn, 3);
    newIiy    =   HVDC.newv(IsOn, 4);          oldIiy    =   HVDC.state(IsOn, 4);
    newIdr    =   HVDC.newv(IsOn, 5);          oldIdr    =   HVDC.state(IsOn, 5);
    newIdi    =   HVDC.newv(IsOn, 6);          oldIdi    =   HVDC.state(IsOn, 6);
    newUc     =   HVDC.newv(IsOn, 7);          oldUc     =   HVDC.state(IsOn, 7);
    newUdro   =   HVDC.newv(IsOn, 8);          oldUdro   =   HVDC.state(IsOn, 8);
    newKar    =   HVDC.newv(IsOn, 9);          oldKar    =   HVDC.state(IsOn, 9);
    newUdr    =   HVDC.newv(IsOn,10);          oldUdr    =   HVDC.state(IsOn,10);
    newmiur   =   HVDC.newv(IsOn,11);          oldmiur   =   HVDC.state(IsOn,11);
    newtanfir =   HVDC.newv(IsOn,12);          oldtanfir =   HVDC.state(IsOn,12);
    newgammar =   HVDC.newv(IsOn,13);          oldgammar =   HVDC.state(IsOn,13);
    newUdio   =   HVDC.newv(IsOn,14);          oldUdio   =   HVDC.state(IsOn,14);
    newKai    =   HVDC.newv(IsOn,15);          oldKai    =   HVDC.state(IsOn,15);
    newUdi    =   HVDC.newv(IsOn,16);          oldUdi    =   HVDC.state(IsOn,16);
    newmiui   =   HVDC.newv(IsOn,17);          oldmiui   =   HVDC.state(IsOn,17);
    newtanfii =   HVDC.newv(IsOn,18);          oldtanfii =   HVDC.state(IsOn,18);
    newgammai =   HVDC.newv(IsOn,19);          oldgammai =   HVDC.state(IsOn,19);
    newalphar =   HVDC.newv(IsOn,20);          oldalphar =   HVDC.state(IsOn,20);
    newalphai =   HVDC.newv(IsOn,21);          oldalphai =   HVDC.state(IsOn,21);
    
    Urx       = HVDC.Urx(IsOn,1);           Ury       = HVDC.Ury(IsOn,1);             Ur = sqrt(Urx.*Urx+Ury.*Ury);
    Uix       = HVDC.Uix(IsOn,1);           Uiy       = HVDC.Uiy(IsOn,1);             Ui = sqrt(Uix.*Uix+Uiy.*Uiy);
    sp        = HVDC.sp(IsOn,1);
    
    HVDC.box(sp+ 1,1) = Urx.*newIrx + Ury.*newIry + newUdr.*newIdr;
    HVDC.box(sp+ 2,1) = Ury.*newIrx - Urx.*newIry + newUdr.*newIdr.*newtanfir;
    HVDC.box(sp+ 3,1) = Uix.*newIix + Uiy.*newIiy - newUdi.*newIdi;
    HVDC.box(sp+ 4,1) = Uiy.*newIix - Uix.*newIiy + newUdi.*newIdi.*newtanfii;
    HVDC.box(sp+ 5,1) = 2*(newUdr+oldUdr) - 2*(newUc+oldUc) - R.*(newIdr+oldIdr) - 2*L.*(newIdr-oldIdr)/hh;
    HVDC.box(sp+ 6,1) = 2*(newUdi+oldUdi) - 2*(newUc+oldUc) + R.*(newIdi+oldIdi) + 2*L.*(newIdi-oldIdi)/hh;
    HVDC.box(sp+ 7,1) = (newIdr+oldIdr) - (newIdi+oldIdi) - 2*C.*(newUc-oldUc)/hh;
    HVDC.box(sp+ 8,1) = newUdro - 3*sqrt(2)/pi*Nbr.*newKar.*Ur;
    HVDC.box(sp+ 9,1) = newKar - Kars;
    HVDC.box(sp+10,1) = newUdr - newUdro.*cos(newalphar) + Nbr.^2.*(3/pi*Xcr+2*Rcr).*newIdr + 2*Nbr.*Vdropr;
    HVDC.box(sp+11,1) = newUdro.*cos(newalphar+newmiur) - newUdro.*cos(newalphar) + 6/pi*Nbr.^2.*Xcr.*newIdr;
    HVDC.box(sp+12,1) = newtanfir.*(cos(2*newalphar)-cos(2*newalphar+2*newmiur)) - (2*newmiur+sin(2*newalphar)-sin(2*newalphar+2*newmiur));
    HVDC.box(sp+13,1) = newalphar + newmiur + newgammar - pi;
    HVDC.box(sp+14,1) = newUdio - 3*sqrt(2)/pi*Nbi.*newKai.*Ui;
    HVDC.box(sp+15,1) = newKai - Kais;
    HVDC.box(sp+16,1) = newUdi + newUdio.*cos(newalphai) - Nbi.^2.*(3/pi*Xci+2*Rci).*newIdi - 2*Nbi.*Vdropi;
    HVDC.box(sp+17,1) = newUdio.*cos(newalphai+newmiui) - newUdio.*cos(newalphai) + 6/pi*Nbi.^2.*Xci.*newIdi;
    HVDC.box(sp+18,1) = newtanfii.*(cos(2*newgammai)-cos(2*newgammai+2*newmiui)) - (2*newmiui+sin(2*newgammai)-sin(2*newgammai+2*newmiui));
    HVDC.box(sp+19,1) = newalphai + newmiui + newgammai - pi;
    HVDC.box(sp+20,1) = 2*Kpr.*(newIdr-oldIdr)/hh - Kir.*(2*Iord-newIdr-oldIdr) - 2*(newalphar-oldalphar)/hh;
    HVDC.box(sp+21,1) = 2*Kpi.*(newgammai-oldgammai)/hh - Kii.*(2*HVDC.gammais(IsOn,1)-newgammai-oldgammai) - 2*(newalphai-oldalphai)/hh;
    
    HVDC.g4n(sp+ 1,1) = newIrx;                                  HVDC.g4n(sp+ 1,2) = newIry;
    HVDC.g4n(sp+ 2,1) = -newIry;                                 HVDC.g4n(sp+ 2,2) = newIrx;
    HVDC.g4n(sp+ 3,3) = newIix;                                  HVDC.g4n(sp+ 3,4) = newIiy;
    HVDC.g4n(sp+ 4,3) = -newIiy;                                 HVDC.g4n(sp+ 4,4) = newIix;
    HVDC.g4n(sp+ 8,1) = -3*sqrt(2)/pi.*newKar.*Nbr.*Urx./Ur;     HVDC.g4n(sp+ 8,2) = -3*sqrt(2)/pi.*newKar.*Nbr.*Ury./Ur;
    HVDC.g4n(sp+14,3) = -3*sqrt(2)/pi.*newKai.*Nbi.*Uix./Ui;     HVDC.g4n(sp+14,4) = -3*sqrt(2)/pi.*newKai.*Nbi.*Uiy./Ui;
    
    y11   = Urx;                                                           I11   = sp+ 1;     J11   = sp+ 1;
    y12   = Ury;                                                           I12   = sp+ 1;     J12   = sp+ 2;
    y15   = newUdr;                                                        I15   = sp+ 1;     J15   = sp+ 5;
    y0110 = newIdr;                                                        I0110 = sp+ 1;     J0110 = sp+10;
    y21   = Ury;                                                           I21   = sp+ 2;     J21   = sp+ 1;
    y22   = -Urx;                                                          I22   = sp+ 2;     J22   = sp+ 2;
    y25   = newUdr.*newtanfir;                                             I25   = sp+ 2;     J25   = sp+ 5;
    y0210 = newIdr.*newtanfir;                                             I0210 = sp+ 2;     J0210 = sp+10;
    y0212 = newUdr.*newIdr;                                                I0212 = sp+ 2;     J0212 = sp+12;
    y33   = Uix;                                                           I33   = sp+ 3;     J33   = sp+ 3;
    y34   = Uiy;                                                           I34   = sp+ 3;     J34   = sp+ 4;
    y36   = -newUdi;                                                       I36   = sp+ 3;     J36   = sp+ 6;
    y0316 = -newIdi;                                                       I0316 = sp+ 3;     J0316 = sp+16;
    y43   = Uiy;                                                           I43   = sp+ 4;     J43   = sp+ 3;
    y44   = -Uix;                                                          I44   = sp+ 4;     J44   = sp+ 4;
    y46   = newUdi.*newtanfii;                                             I46   = sp+ 4;     J46   = sp+ 6;
    y0416 = newIdi.*newtanfii;                                             I0416 = sp+ 4;     J0416 = sp+16;
    y0418 = newUdi.*newIdi;                                                I0418 = sp+ 4;     J0418 = sp+18;
    y55   = -R-2*L/hh;                                                     I55   = sp+ 5;     J55   = sp+ 5;
    y57   = -2*unit;                                                       I57   = sp+ 5;     J57   = sp+ 7;
    y0510 = 2*unit;                                                        I0510 = sp+ 5;     J0510 = sp+10;
    y66   = R+2*L/hh;                                                      I66   = sp+ 6;     J66   = sp+ 6;
    y67   = -2*unit;                                                       I67   = sp+ 6;     J67   = sp+ 7;
    y0616 = 2*unit;                                                        I0616 = sp+ 6;     J0616 = sp+16;
    y75   = unit;                                                          I75   = sp+ 7;     J75   = sp+ 5;
    y76   = -unit;                                                         I76   = sp+ 7;     J76   = sp+ 6;
    y77   = -2*C/hh;                                                       I77   = sp+ 7;     J77   = sp+ 7;
    y88   = unit;                                                          I88   = sp+ 8;     J88   = sp+ 8;
    y89   = -3*sqrt(2)/pi*Nbr.*Ur;                                         I89   = sp+ 8;     J89   = sp+ 9;
    y99   = unit;                                                          I99   = sp+ 9;     J99   = sp+ 9;
    y1005 = Nbr.^2.*(3/pi*Xcr+2*Rcr);                                      I1005 = sp+10;     J1005 = sp+ 5;
    y1008 = -cos(newalphar);                                               I1008 = sp+10;     J1008 = sp+ 8;
    y1010 = unit;                                                          I1010 = sp+10;     J1010 = sp+10;
    y1020 = newUdro.*sin(newalphar);                                       I1020 = sp+10;     J1020 = sp+20;
    y1105 = 6/pi*Nbr.^2.*Xcr;                                              I1105 = sp+11;     J1105 = sp+ 5;
    y1108 = cos(newalphar+newmiur)-cos(newalphar);                         I1108 = sp+11;     J1108 = sp+ 8;
    y1111 = -newUdro.*sin(newalphar+newmiur);                              I1111 = sp+11;     J1111 = sp+11;
    y1120 = -newUdro.*sin(newalphar+newmiur)+newUdro.*sin(newalphar);      I1120 = sp+11;     J1120 = sp+20;
    y1211 = 2*newtanfir.*sin(2*newalphar+2*newmiur)-2*unit+2*cos(2*newalphar+2*newmiur);
                                                                           I1211 = sp+12;     J1211 = sp+11;
    y1212 = cos(2*newalphar)-cos(2*newalphar+2*newmiur);                   I1212 = sp+12;     J1212 = sp+12;
    y1220 = 2*newtanfir.*(-sin(2*newalphar)+sin(2*newalphar+2*newmiur))-2*cos(2*newalphar)+2*cos(2*newalphar+2*newmiur);
                                                                           I1220 = sp+12;     J1220 = sp+20;
    y1311 = unit;                                                          I1311 = sp+13;     J1311 = sp+11;
    y1313 = unit;                                                          I1313 = sp+13;     J1313 = sp+13;
    y1320 = unit;                                                          I1320 = sp+13;     J1320 = sp+20;
    y1414 = unit;                                                          I1414 = sp+14;     J1414 = sp+14;
    y1415 = -3*sqrt(2)/pi*Nbi.*Ui;                                         I1415 = sp+14;     J1415 = sp+15;
    y1515 = unit;                                                          I1515 = sp+15;     J1515 = sp+15;
    y1606 = -Nbi.^2.*(3/pi*Xci+2*Rci);                                     I1606 = sp+16;     J1606 = sp+ 6;
    y1614 = cos(newalphai);                                                I1614 = sp+16;     J1614 = sp+14;
    y1616 = unit;                                                          I1616 = sp+16;     J1616 = sp+16;
    y1621 = -newUdio.*sin(newalphai);                                      I1621 = sp+16;     J1621 = sp+21;
    y1706 = 6/pi*Nbi.^2.*Xci;                                              I1706 = sp+17;     J1706 = sp+ 6;
    y1714 = cos(newalphai+newmiui)-cos(newalphai);                         I1714 = sp+17;     J1714 = sp+14;
    y1717 = -newUdio.*sin(newalphai+newmiui);                              I1717 = sp+17;     J1717 = sp+17;
    y1721 = -newUdio.*sin(newalphai+newmiui)+newUdio.*sin(newalphai);      I1721 = sp+17;     J1721 = sp+21;
    y1817 = 2*newtanfii.*sin(2*newgammai+2*newmiui)-2*unit+2*cos(2*newgammai+2*newmiui); 
                                                                           I1817 = sp+18;     J1817 = sp+17;
    y1818 = cos(2*newgammai)-cos(2*newgammai+2*newmiui);                   I1818 = sp+18;     J1818 = sp+18;
    y1819 = 2*newtanfii.*(-sin(2*newgammai)+sin(2*newgammai+2*newmiui))-2*cos(2*newgammai)+2*cos(2*newgammai+2*newmiui);
                                                                           I1819 = sp+18;     J1819 = sp+19;
    y1917 = unit;                                                          I1917 = sp+19;     J1917 = sp+17;
    y1919 = unit;                                                          I1919 = sp+19;     J1919 = sp+19;
    y1921 = unit;                                                          I1921 = sp+19;     J1921 = sp+21;
    y2020 = -2/hh*unit;                                                    I2020 = sp+20;     J2020 = sp+20;
    y20x1 = 2*Kpr/hh+Kir;                                                  I20x1 = sp+20;     J20x1 = sp+ 5;
    y20x2 = [];                                                            I20x2 = [];        J20x2 = [];
    y2121 = -2/hh*unit;                                                    I2121 = sp+21;     J2121 = sp+21;
    y21x1 = 2*Kpi/hh+Kii;                                                  I21x1 = sp+21;     J21x1 = sp+19;
    y21x2 = [];                                                            I21x2 = [];        J21x2 = [];

%测试南方电网
%     HVDC.box(HVDC.sp+20,1) = newalphar - HVDC.alphar;
%     HVDC.box(HVDC.sp+21,1) = newalphai - HVDC.alphai;
%     y2020 = unit;                                                          I2020 = sp+20;     J2020 = sp+20;
%     y20x1 = [];                                                            I20x1 = [];             J20x1 = [];
%     y20x2 = [];                                                            I20x2 = [];             J20x2 = [];
%     y2121 = unit;                                                          I2121 = sp+21;     J2121 = sp+21;
%     y21x1 = [];                                                            I21x1 = [];             J21x1 = [];
%     y21x2 = [];                                                            I21x2 = [];             J21x2 = [];

    %% 故障瞬间，直流系统的被控变量为alphar和alphai，它们不发生变化
    if faultflag.postchange==1
        HVDC.box(sp+ 5,1) = 0;
        y55   = unit;                                                      I55   = sp+ 5;     J55   = sp+ 5;
        y57   = [];                                                        I57   = [];        J57   = [];
        y0510 = [];                                                        I0510 = [];        J0510 = [];
        HVDC.box(sp+ 6,1) = 0;
        y66   = unit;                                                      I66   = sp+ 6;     J66   = sp+ 6;
        y67   = [];                                                        I67   = [];        J67   = [];
        y0616 = [];                                                        I0616 = [];        J0616 = [];
        HVDC.box(sp+ 7,1) = 0;
        y75   = [];                                                        I75   = [];        J75   = [];
        y76   = [];                                                        I76   = [];        J76   = [];
        y77   = unit;                                                      I77   = sp+ 7;     J77   = sp+ 7;
        
        HVDC.box(sp+11,1) = 0;
        y1105 = [];                                                        I1105 = [];        J1105 = [];
        y1108 = [];                                                        I1108 = [];        J1108 = [];
        y1111 = unit;                                                      I1111 = sp+11;     J1111 = sp+11;
        y1120 = [];                                                        I1120 = [];        J1120 = [];
%         HVDC.box(sp+12,1) = 0;
%         y1211 = [];                                                        I1211 = [];        J1211 = [];
%         y1212 = unit;                                                      I1212 = sp+12;     J1212 = sp+12;
%         y1220 = [];                                                        I1220 = [];        J1220 = [];
        HVDC.box(sp+13,1) = 0;
        y1311 = [];                                                        I1311 = [];        J1311 = [];
        y1313 = unit;                                                      I1313 = sp+13;     J1313 = sp+13;
        y1320 = [];                                                        I1320 = [];        J1320 = [];
        HVDC.box(sp+17,1) = 0;
        y1706 = [];                                                        I1706 = [];        J1706 = [];
        y1714 = [];                                                        I1714 = [];        J1714 = [];
        y1717 = unit;                                                      I1717 = sp+17;     J1717 = sp+17;
        y1721 = [];                                                        I1721 = [];        J1721 = [];
%         HVDC.box(sp+18,1) = 0;
%         y1817 = [];                                                        I1817 = [];        J1817 = [];
%         y1818 = unit;                                                      I1818 = sp+18;     J1818 = sp+18;
%         y1819 = [];                                                        I1819 = [];        J1819 = [];
        HVDC.box(sp+19,1) = 0;
        y1917 = [];                                                        I1917 = [];        J1917 = [];
        y1919 = unit;                                                      I1919 = sp+19;     J1919 = sp+19;
        y1921 = [];                                                        I1921 = [];        J1921 = [];

        HVDC.box(sp+20,1) = 0;
        y2020 = unit;                                                      I2020 = sp+20;     J2020 = sp+20;
        y20x1 = [];                                                        I20x1 = [];        J20x1 = [];
        y20x2 = [];                                                        I20x2 = [];        J20x2 = [];
        HVDC.box(sp+21,1) = 0;
        y2121 = unit;                                                      I2121 = sp+21;     J2121 = sp+21;
        y21x1 = [];                                                        I21x1 = [];        J21x1 = [];
        y21x2 = [];                                                        I21x2 = [];        J21x2 = [];
        fprintf('换路\n');
    end
    
    %% 整流侧逆变侧角度越限处理
%     if ((faultflag.postchange==0)&&(~isempty(IsLimited.Ar)))
%         no = IsLimited.Ar;                                sp = HVDC.sp(no);
%         HVDC.box(sp+20,1) = newalphar(no,1) - HVDC.Aminr(no,1);
%         y2020(no,1) = unit(no,1);
%         y20x1(no,1) = 0;
%     end
%     if ((faultflag.postchange==0)&&(~isempty(IsLimited.Gi)))
%         no = IsLimited.Gi;                                sp = HVDC.sp(no);
%         HVDC.box(sp+21,1) = newgammai(no) - HVDC.gammais(no);
%         y2121(no,1) = 0;
%         y21x1(no,1) = unit(no,1);
%     end
    
    if (faultflag.postchange==0)
        if ~isempty(IsLimited.Ar)
            index = intersect(IsOn,IsLimited.Ar);
            no = IsOn(index,1);                           sp = HVDC.sp(no);
            HVDC.box(sp+20,1) = newalphar(index,1) - HVDC.Aminr(no,1);
            y2020(index,1) = unit(index,1);
            y20x1(index,1) = 0;
        end
        if ~isempty(IsLimited.Gi)
            index = intersect(IsOn,IsLimited.Gi);
            no = IsOn(index,1);                           sp = HVDC.sp(no);
            HVDC.box(sp+21,1) = newgammai(index,1) - HVDC.Gimin(no,1);
            y2121(index,1) = 0;
            y21x1(index,1) = unit(index,1);
        end
    end
    
    
    %% 闭锁直流线路雅克比矩阵处理
    IsBlock = find(HVDC.IsBlock == 1);
    if ~isempty(IsBlock)
        n = [];      nn = [];
        for i = 1:size(IsBlock,1)
            sp  = HVDC.sp(IsBlock(i,1),1);
            tsn = HVDC.tsn(IsBlock(i,1),1);
            nn(:,1) = [sp:sp+tsn]';
            n  = [n;nn];
        end
        unit = ones(size(n,1),1);
    else
        n    = [];
        unit = [];
    end
    %% 构造直流系统的雅克比矩阵
    dcyI = cat(1,dcyI,I11,I12,I15,I0110,I21,I22,I25,I0210,I0212,I33,I34,I36,I0316,I43,I44,I46,I0416,I0418,I55,I57,I0510,I66,I67,I0616,...
                      I75,I76,I77,I88,I89,I99,I1005,I1008,I1010,I1020,I1105,I1108,I1111,I1120,I1211,I1212,I1220,I1311,I1313,I1320,...
                      I1414,I1415,I1515,I1606,I1614,I1616,I1621,I1706,I1714,I1717,I1721,I1817,I1818,I1819,I1917,I1919,I1921,...
                      I2020,I20x1,I20x2,I2121,I21x1,I21x2,n);
    dcyJ = cat(1,dcyJ,J11,J12,J15,J0110,J21,J22,J25,J0210,J0212,J33,J34,J36,J0316,J43,J44,J46,J0416,J0418,J55,J57,J0510,J66,J67,J0616,...
                      J75,J76,J77,J88,J89,J99,J1005,J1008,J1010,J1020,J1105,J1108,J1111,J1120,J1211,J1212,J1220,J1311,J1313,J1320,...
                      J1414,J1415,J1515,J1606,J1614,J1616,J1621,J1706,J1714,J1717,J1721,J1817,J1818,J1819,J1917,J1919,J1921,...
                      J2020,J20x1,J20x2,J2121,J21x1,J21x2,n);
    dcyV = cat(1,dcyV,y11,y12,y15,y0110,y21,y22,y25,y0210,y0212,y33,y34,y36,y0316,y43,y44,y46,y0416,y0418,y55,y57,y0510,y66,y67,y0616,...
                      y75,y76,y77,y88,y89,y99,y1005,y1008,y1010,y1020,y1105,y1108,y1111,y1120,y1211,y1212,y1220,y1311,y1313,y1320,...
                      y1414,y1415,y1515,y1606,y1614,y1616,y1621,y1706,y1714,y1717,y1721,y1817,y1818,y1819,y1917,y1919,y1921,...
                      y2020,y20x1,y20x2,y2121,y21x1,y21x2,unit);
    %% 矩阵预算
    HVDC.box = -HVDC.box;
    dcy = sparse(dcyI, dcyJ, dcyV);
    dim = size(dcy,1);

    inv_dcy = dcy\speye(dim);
    HVDC.g4n = inv_dcy*HVDC.g4n;
    HVDC.box = inv_dcy*HVDC.box;
    
    HVDC.Irx = HVDC.newv(:,1);                   HVDC.Iry = HVDC.newv(:,2); %整流器，直流系统注入交流系统电流
    HVDC.Iix = HVDC.newv(:,3);                   HVDC.Iiy = HVDC.newv(:,4); %逆变器，直流系统注入交流系统电流
end

return;