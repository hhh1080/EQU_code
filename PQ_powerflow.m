function [bus, iter]=PQ_powerflow(bus,branch,Y)
threshold=1e-10;
j=sqrt(-1);
S = (bus(:,4) - bus(:,6))+j*(bus(:,5)- bus(:,7));
pg = real(S);                                      qg = imag(S);
bus_num = size(bus,1);
v = bus(:,2);                                      o = 0*pi/180*ones(bus_num,1);
flags=0;                       
refbus	= find(bus(:,12) == 3 );                   pvbus	= find(bus(:,12) == 2 );
nl = size(branch, 1);                              nb = size(bus, 1);  
f = branch(:, 1);	                               t = branch(:, 2);							
Cf = sparse(f, 1:nl, ones(nl, 1), nb, nl);         Ct = sparse(t, 1:nl, ones(nl, 1), nb, nl);
gl =find(f==t);            Cf(:,gl)=0;             Ct(:,gl)=0;
X   = branch(:,4);
X(X==0)=1e-10;
Bs   = 1./X;
Btt = Bs;                                          Bff = Bs;
Bft = - Bs;                                        Btf = - Bs;
B= Cf * spdiags(Bff, 0, nl, nl) * Cf' + ...
    Cf * spdiags(Bft, 0, nl, nl) * Ct' + ...
    Ct * spdiags(Btf, 0, nl, nl) * Cf' + ...
    Ct * spdiags(Btt, 0, nl, nl) * Ct';
BB=-imag(Y);
while flags==0
    iter=0;
    Jp=B;                Jq=BB;
    Jp(refbus,:)=0;      Jq(refbus,:)=0;            Jq(pvbus,:)=0;
    dJvv =sparse(refbus,refbus,1,nb,nb);            dJpv =sparse(pvbus,pvbus,1,nb,nb);
    Jp  = Jp  + dJvv;                               Jq =  Jq + dJvv +dJpv;
    while 1
        U=v.*(cos(o)+j*sin(o));                     Pe=real(U.*conj(Y*U));
        dP=(pg-Pe)./v;                              dP(refbus)=0;
        do=Jp\dP;
        o=o+do;
        if (any(abs(dP)>threshold))
            kp=0;
        else
            kp=1;
            if kq==1 
                break; 
            end
        end
        U=v.*(cos(o)+j*sin(o));                     Qe=imag(U.*conj(Y*U));
        dQ=(qg-Qe)./v;         dQ(refbus)=0;        dQ(pvbus)=0;
        dv=Jq\dQ;
        v=v+dv;
        if (any(abs(dQ)>threshold))
            kq=0;
        else
            kq=1;
            if kp==1 
                break; 
            end
        end    
        iter=iter+1;  
        if iter>100
            error('the load flow does not converge.');
        end
    end
    U=v.*(cos(o)+j*sin(o));                        Qe=imag(U.*conj(Y*U));
%    gq   =Qe+bus(:,7);                             gqmax=bus(:,8);
     flags=1;
%     index1=find(gq>gqmax);                         index2=find(bus(:,12)==2);          index=intersect(index1,index2);
 %   if size(index,1)==0
  %      flags=1;
   % else
    %%   bus(index,12)=1;                    pvbus	= find(bus(:,12) == 2 );
      %  bus(index,5) = bus(index,8);        qg = bus(:,5) - bus(:,7);
       % kp=0;                               kq=0;
%    end
end
bus(:,2)  = v;           bus(:,3) = o;    U=v.*(cos(o)+j*sin(o));
Se=U.*conj(Y*U);         Pe=real(Se);     Qe=imag(Se);
bus(:,15)    = Pe;                        bus(:,16)    = Qe;
bus(:,17) = Pe + j.*Qe + bus(:,6) + j.*bus(:,7);