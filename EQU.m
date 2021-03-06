function [EQUbus,EQUbranch,EQUgen,EQUY,EQUfault] = EQU(bus,branch,Y,gen,savebus,savebranch,group,fault)

busnum=size(bus,1);
%找到与group关联的节点
groupnum=size(group,1);%%groupnum为机群数
[groupnum,groupgennum]=size(group);
busno=bus(:,1);
mbusno=busno;
EQUY=Y;
V=bus(:,2).*cos(bus(:,3))+sqrt(-1)*bus(:,2).*sin(bus(:,3));
mV=V;
volt=bus(:,14);

EQUbus=bus;
EQUbranch=[];
tempbranch=branch;
EQUgen=gen;
%% 等值发电机母线化简
for ii=1:groupnum
    equbusno=busnum+ii;  %等值母线编号
   %找到要消去节点编号
    CBUS=group(ii,find(group(ii,2:groupgennum))+1);
    CBUS=CBUS';
   %计算等值母线电压  
    vemagn=0;
    veangle=0;
    for i=1:size(CBUS,1)
        vemagn=vemagn+bus(CBUS(i,1),2);
        veangle=veangle+bus(CBUS(i,1),3);
    end
    vemagn=vemagn/size(CBUS,1);%%用电压的平均值来表示等值母线电压
    veangle=veangle/size(CBUS,1);%%用角度的平均值来表示等值母线角度
    ve=vemagn*cos(veangle)+sqrt(-1)*vemagn*sin(veangle);
    
    %找到与相应等值集群相关联的节点群（就是与需要等值的节点有线路连接的节点）
    RBUS=[];
    for jj=2:groupgennum
        if group(ii,jj)~=0
            tempR=find(Y(group(ii,jj),:));
            tempR=setdiff(tempR,group(ii,jj)); 
            RBUS=union(RBUS,tempR);
        end
    end
    RBUS=RBUS';
    %寻找完毕，存在向量RBUS中  
    
    %求发电机母线简化后的系统导纳矩阵
    VC=V(CBUS);
    VR=V(RBUS);
    YRC=Y(RBUS,CBUS); 
   % YRE=YRC*(VC/ve);
    YRE=YRC*abs(VC/ve);    %消去移相器
    YCR=Y(CBUS,RBUS);
   % YER=(VC/ve)'*YCR; % '对于实数矩阵是求转置；对于复数矩阵是求共轭转置
    YER=(abs(VC/ve))'*YCR; %消去移相器
    YCC=Y(CBUS,CBUS);
    YEE=(VC/ve)'*YCC*(VC/ve);    
    YRR=diag(Y(RBUS,RBUS));  %n*1列向量
%    YRR_new=YRR+YRC*ones(size(CBUS,1),1)-YRE; %更新对角阵,这样更新后不正确。
    YRR_new=YRR;
    %消去移相器，特殊处理
    yR0=(ve./VR).*(YRC*(abs(VC/ve).*(cos(atan2(imag(VC/ve),real(VC/ve)))+sqrt(-1)*sin(atan2(imag(VC/ve),real(VC/ve)))-1)));
    yE0=(conj(VR/ve))'*(YRC*(abs(VC/ve).*(cos(atan2(imag(VC/ve),real(VC/ve)))-sqrt(-1)*sin(atan2(imag(VC/ve),real(VC/ve)))-1)));
    YRR_new=YRR_new+yR0;
    YEE=YEE+yE0;
    
    %更新导纳矩阵
    [mbusno,mI]=setdiff(mbusno,CBUS); %去除掉要消去的节点
    mV=mV(mI);  %电压值中也取消要消去的节点电压
    EQUY=EQUY(mI,mI); %要保留的导纳矩阵    
    EQUbus=EQUbus(mI,1:15);  %保留节点的bus参数-----new    
    dim_I=size(mI,1);  %目前节点个数
    for jj=1:size(RBUS,1)
        R=find(mbusno==RBUS(jj));
        EQUY(R,dim_I+1)=YRE(jj,1);
        EQUY(dim_I+1,R)=YER(1,jj);
        EQUY(R,R)=YRR_new(jj,1);
    end
    EQUY(dim_I+1,dim_I+1)=YEE;
    mbusno(dim_I+1)=equbusno; %增加化简节点
   
    %更新EQUbus
    newbus=zeros(1,15);
    newbus(1,1)=equbusno;newbus(1,2)=abs(ve); newbus(1,3)=atan2(imag(ve),real(ve));newbus(1,4)=sum(bus(CBUS,4));
    newbus(1,5)=sum(bus(CBUS,5)); newbus(1,6)=sum(bus(CBUS,6)); newbus(1,7)=sum(bus(CBUS,7)); newbus(1,8)=sum(bus(CBUS,8)); 
    newbus(1,9)=sum(bus(CBUS,9)); newbus(1,10)=sum(bus(CBUS,10)); newbus(1,11)=sum(bus(CBUS,11)); newbus(1,12)=sum(bus(CBUS,12)); 
    
    if ~isempty(find(bus(CBUS,13)==3))
     %   msgbox('当前发电机组中有平衡机，将等值机设为平衡机','提醒');
        newbus(1,13)=3; 
    elseif ~isempty(find(bus(CBUS,13)==2))
     %   msgbox('当前发电机组中存在PV节点，等值机类型为PV节点','提醒');
        newbus(1,13)=2; 
    else
        newbus(1,13)=1;
    end
    newbus(1,14)=(bus(CBUS(1,1),14));
    EQUbus=[EQUbus;newbus];  %增加的等值节点加入
    
    mV(dim_I+1)=ve; %增加化简节点电压
    
    %发电机等值母线化简后的等值母线是要保留的节点
    savebus=[savebus;equbusno];
end
% 等值发电机母线化简结束
%% 网络化简
% 辨别保留节点、边界节点和消去节点(mbusno实际代表母线名)
% mbusno为同调等值后的全部节点号，savebus是保留母线的节点号
[throwbus,throwno]=setdiff(mbusno,savebus); 
[savebus,saveno]=setdiff(mbusno,throwbus);

thrownum=size(throwbus,1); %消去节点个数
borderbus=[];
for ii=1:thrownum
    tempbus=mbusno(saveno(find(EQUY(throwno(ii),saveno)~=0)));
    tempbus=setdiff(tempbus,mbusno(throwno));
    borderbus=union(borderbus,tempbus);
end

%将borderbus从savebus中去掉
%borderno,saveno,throwno都是在mbusno中的索引。
[innerbus,innerno]=setdiff(savebus,borderbus);
innerno=saveno(innerno);
[borderbus,borderno]=setdiff(savebus,innerbus);
borderno=saveno(borderno);

%求电流
I=EQUY*mV;

innerV=mV(innerno);
borderV=mV(borderno);
throwV=mV(throwno);

YII=EQUY(innerno,innerno);
YIB=EQUY(innerno,borderno);
YBI=EQUY(borderno,innerno);
YBB=EQUY(borderno,borderno);
YBD=EQUY(borderno,throwno);
YDB=EQUY(throwno,borderno);
YDD=EQUY(throwno,throwno);
YDI=EQUY(throwno,innerno);

YBB_new=YBB-YBD*inv(YDD)*YDB;
%更新EQUY、EQUV
EQUY=[YII YIB;YBI YBB_new];
EQUV=[innerV;borderV];
%消去节点转移到边界节点的等值负荷
IB_new=I(borderno)-YBD*inv(YDD)*I(throwno);
I_new=[I(innerno);IB_new];
s_new=[innerV.*conj(I(innerno));borderV.*conj(IB_new)];
load_new=borderV.*conj(YBD*inv(YDD)*I(throwno));
EQUbus(borderno,6)=EQUbus(borderno,6)+real(load_new);
EQUbus(borderno,7)=EQUbus(borderno,7)+imag(load_new);
equbus=[mbusno(innerno);mbusno(borderno)];
EQUbus=[EQUbus(innerno,:);EQUbus(borderno,:)];

%% 输出潮流数据
equbusnum=size(equbus,1);
EQUbranch=savebranch;
tempEQUY=EQUY;
savebranchnum=size(savebranch,1);

for ii=1:savebranchnum
    ibus=savebranch(ii,1);
    jbus=savebranch(ii,2);
    ino=find(equbus==ibus);
    jno=find(equbus==jbus);
    if savebranch(ii,8)==1
        tempEQUY(ino,jno)=tempEQUY(ino,jno)+1/(savebranch(ii,3)+sqrt(-1)*savebranch(ii,4));
        tempEQUY(jno,ino)=tempEQUY(jno,ino)+1/(savebranch(ii,3)+sqrt(-1)*savebranch(ii,4));
        tempEQUY(ino,ino)=tempEQUY(ino,ino)-sqrt(-1)*savebranch(ii,5)-1/(savebranch(ii,3)+sqrt(-1)*savebranch(ii,4));
        tempEQUY(jno,jno)=tempEQUY(jno,jno)-sqrt(-1)*savebranch(ii,6)-1/(savebranch(ii,3)+sqrt(-1)*savebranch(ii,4));
    elseif savebranch(ii,8)==2
        tempEQUY(ino,jno)=tempEQUY(ino,jno)+1/((savebranch(ii,3)+sqrt(-1)*savebranch(ii,4))*savebranch(ii,7));
        tempEQUY(jno,ino)=tempEQUY(jno,ino)+1/((savebranch(ii,3)+sqrt(-1)*savebranch(ii,4))*savebranch(ii,7));
        tempEQUY(ino,ino)=tempEQUY(ino,ino)-1/((savebranch(ii,3)+sqrt(-1)*savebranch(ii,4))*savebranch(ii,7)*savebranch(ii,7));
        tempEQUY(jno,jno)=tempEQUY(jno,jno)-1/(savebranch(ii,3)+sqrt(-1)*savebranch(ii,4));
    end
end

diagY=diag(tempEQUY);
for ii=1:equbusnum
    nonzeroelementno=find(abs(tempEQUY(ii,:))>1e-5);
    nonzeroelementno=setdiff(nonzeroelementno,ii);  %与ii相关联的节点编号
    if ~isempty(nonzeroelementno)
        for jj=1:size(nonzeroelementno',1)  %在用到size的时候，一定要注意是行向量还是列向量。
            diagY(ii)=diagY(ii)+tempEQUY(ii,nonzeroelementno(jj));
        end  
    end
end


for ii=1:equbusnum
    nonzeroelementno=find(abs(tempEQUY(ii,:))>1e-6);
    for jj=1:size(nonzeroelementno',1)
        if nonzeroelementno(jj)>ii
            addline=zeros(1,8);
            zs=-1/tempEQUY(ii,nonzeroelementno(jj));
            addline(1,1)=equbus(ii,1);
            addline(1,2)=equbus(nonzeroelementno(jj),1);
            addline(1,3)=real(zs);
            addline(1,4)=imag(zs);
            addline(1,5)=0;
            addline(1,6)=0;
            addline(1,7)=1;
            addline(1,8)=1;
            EQUbranch=[EQUbranch;addline];
        end
    end
end
for ii=1:size(diagY,1)
    if abs(diagY(ii,1))<1e-8
        diagY(ii,1)=0;
    end
end
EQUbus(:,11)=EQUbus(:,11)+real(diagY);
EQUbus(:,12)=EQUbus(:,12)+imag(diagY);
%% 节点重新排序
EQUbranchnum=size(EQUbranch,1);
for ii=1:EQUbranchnum
   EQUbranch(ii,1)=find(equbus==EQUbranch(ii,1));
   EQUbranch(ii,2)=find(equbus==EQUbranch(ii,2));
end
for ii=1:size(EQUbus,1)
    EQUbus(ii,1)=ii;
end

%% 发电机聚合 目前只考虑经典模型，时间常数也需要修改，否则误差很大呢，初始角度误差就很大。
for ii=1:size(group,1)  
    [genm,genn]=size(gen);
    tempgen=zeros(1,genn);
    CBUS=group(ii,find(group(ii,2:groupgennum))+1);
    CBUS=CBUS';
    [savegen,savegenno]=setdiff(EQUgen(:,2),CBUS);
    EQUgen=EQUgen(savegenno,:);
    tempgen(1,1)=size(EQUgen,1)+1;
    tempgen(1,2)=busnum+ii;
    for jj=1:size(CBUS,1)
        tempgen(1,3)=tempgen(1,3)+gen(find(gen(:,2)==CBUS(jj,1)),3);
        tempgen(1,5)=tempgen(1,5)+1/(gen(find(gen(:,2)==CBUS(jj,1)),5));
        tempgen(1,6)=tempgen(1,6)+1/(gen(find(gen(:,2)==CBUS(jj,1)),6));
        tempgen(1,7)=tempgen(1,7)+1/(gen(find(gen(:,2)==CBUS(jj,1)),7));
        tempgen(1,8)=tempgen(1,8)+1/(gen(find(gen(:,2)==CBUS(jj,1)),8));  
        tempgen(1,9)=tempgen(1,9)+1/(gen(find(gen(:,2)==CBUS(jj,1)),9)); 
        tempgen(1,10)=tempgen(1,10)+1/(gen(find(gen(:,2)==CBUS(jj,1)),10)); 
    end
    for jj=5:10
        tempgen(1,jj)=1/tempgen(1,jj);
    end
    tempgen(1,15)=1;
    tempgen(1,16)=0;
    EQUgen=[EQUgen;tempgen];
end

for ii=1:size(EQUgen,1)
    EQUgen(ii,2)=find(equbus==EQUgen(ii,2));
    EQUgen(ii,1)=ii;
end
EQUgen(:,15)=1;  %发电机类型先设为1，即经典模型
%% 故障一致
EQUfault=fault;
for ii=1:size(EQUfault,1)
    EQUfault(ii,2)=find(equbus==EQUfault(ii,2));
    EQUfault(ii,3)=find(equbus==EQUfault(ii,3));
end

%% 整理并输出等值结果
a=[1 2 3;4 5 6];
equfid=fopen('caseequ.m','wt');
fprintf(equfid,'%s\n','function [bus, branch, Dbus, Dbranch, OLTC, gen, avr, pss, gov, DT, fault, system]=caseequ');
fprintf(equfid,'%s\n','system.MVABASE = 100;');
fprintf(equfid,'%s\n','system.fhz     = 50; ');
fprintf(equfid,'%s\n','system.Tmax    = 6; ');
fprintf(equfid,'%s\n','system.hh      = 1/system.fhz;');
fprintf(equfid,'%s\n','bus=[');
fprintf(equfid,'%4d\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%4d\t%8f\t%8f\n',EQUbus');
fprintf(equfid,'];\n');
fprintf(equfid,'%s\n','branch=[');
fprintf(equfid,'%8d\t%8d\t%8f\t%8f\t%8f\t%8f\t%8f\t%8d\n',EQUbranch');
fprintf(equfid,'];\n');
fprintf(equfid,'%s\n','Dbus=[];');
fprintf(equfid,'%s\n','Dbranch=[];');
fprintf(equfid,'%s\n','OLTC=[];');
fprintf(equfid,'%s\n','gen=[');
fprintf(equfid,'%8d\t%8d\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8d\t%8f\t%8f\n',EQUgen');
fprintf(equfid,'%s\n','];')
fprintf(equfid,'%s\n','avr=[];');
fprintf(equfid,'%s\n','pss=[];');
fprintf(equfid,'%s\n','gov=[];');
fprintf(equfid,'%s\n','DT=[];');
fprintf(equfid,'%s\n','fault=[');
fprintf(equfid,'%8d\t%8d\t%8d\t%8d\t%8f\t%8f\t%8f\t%8f\t%8f\t%8f\t%8d\n',EQUfault');
fprintf(equfid,'%s\n','];');
fclose(equfid);
msgbox('处理结束','提醒');  


