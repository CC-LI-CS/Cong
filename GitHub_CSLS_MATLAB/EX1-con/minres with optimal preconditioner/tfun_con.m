function y = tfun_con(eigs_P4,tv,N,M,Ms)
%% (1) test inv(kron(F*Da,kron(S,S))).'*ones(Ms*N,1)---------------------------the first step
%% test kron(S,S)*Tv
Tv=tv(1:Ms);
TT=[];
for j=1:N
    x=reshape(Tv(:,j),M,M);
    t1t=sqrt(2/(M+1))*dst(sqrt(2/(M+1))*dst(x).').';
    TT(:,j)=t1t(:);
end
Vin=TT; 

Tv1=tv(Ms+1:end);
TT=[];
for j=1:N
    x=reshape(Tv1(:,j),M,M);
    t1t=sqrt(2/(M+1))*dst(sqrt(2/(M+1))*dst(x).').';
    TT(:,j)=t1t(:);
end
Vin1=TT; 

Vin=[Vin;Vin1];
%J*N
% %% Vin*(F'inv(Da))T=(Vin*F')inv(Da) :J*N N*N
% % Vin*F' ===== (F')T*(Vin)T=F'(Vin)T
% %%%test Vin * inv(Da)
% % tv=reshape(v1,N,M-1);
% Vout=(Vin'.*(1./da))';
% Vin=Vout.';
% TT=[];
% for i=1:Ms
%     TT(:,i)=sqrt(N)*ifft(Vin(:,i));
% end
% Vout=TT.';
% %% checked  
% %% (2) test EIG*inv(kron(F*Da,kron(S,S))).'*ones(Ms*N,1)---------------------------the second step
% v=zeros(N,1);
% w0=zeros(N*Ms,1);
% for i=1:N
%     v(i)=1;
%     lam=fft(da.*v);
% %     lam=sqrt(N)*F*Da*v;
%     tI=cT(i)*ones(Ms,1);
%     y=kron(lam,tI);
%     y=y+w0;
%     w0=y;
%     v=zeros(N,1);
% end
% % size(y)
% lamL1=h^(-2)*eigs';
% lam=kron(ones(N,1),lamL1)+y;
vin=eigs_P4.^(-1).*Vin(:);

%% (3) test  inv(kron(inv(Da)*F',kron(S,S))).'  *   < EIG*inv(kron(F*Da,kron(S,S))).'*ones(Ms*N,1)>---------------------------the third step

%kron(S,S)*vin

Vin=reshape(vin(1:Ms),Ms,N);
TT=[];
for j=1:N
    x=reshape(Vin(:,j),M,M);
    t1t=sqrt(2/(M+1))*dst(sqrt(2/(M+1))*dst(x).').';
    TT(:,j)=t1t(:);
end
Vout=TT; %J*N

Vin1=reshape(vin(Ms+1:end),Ms,N);
TT=[];
for j=1:N
    x=reshape(Vin1(:,j),M,M);
    t1t=sqrt(2/(M+1))*dst(sqrt(2/(M+1))*dst(x).').';
    TT(:,j)=t1t(:);
end
Vout1=TT; %J*N
Vout=[Vout;Vout1];

% %% Vin*(FDa)=(Vin*F)(Da) :J*N N*N
% % Vin*F' ===== (F')T*(Vin)T=F'(Vin)T
% 
% Vin=Vout.';
% TT=[];
% for i=1:Ms
%     TT(:,i)=1./sqrt(N)*fft(Vin(:,i));
% end
% Vin=TT.';
% 
% %%%test Vin * Da
% % tv=reshape(v1,N,M-1);
% Vout=(Vin'.*(da))';
% 
% %% %% finished  Ca'^(-0.5)*ones(Ms*N,1);
% %% next comput Ca^(-0.5)
% %% (4) test  inv(kron(F*Da,kron(S,S))) * Vout  < EIG*inv(kron(F*Da,kron(S,S))).'*ones(Ms*N,1)>---------------------------the third step
% % tx=Ca^(-0.5)*ones(Ms*N,1);
% % testtx=inv(kron(F*Da,kron(S,S)))*diag(lam.^(-0.5))*inv(kron(inv(Da)*F',kron(S,S)))*ones(Ms*N,1);
% % norm(tx-testtx)
% 
% %% (4) compute inv(kron(inv(Da)*F',kron(S,S)))*Vout;
% %% copute kron(F*Da,kron(S,S));
% % %% the fourth step compute kron(F*Da,kron(S,S))*vin
% %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^6
% % vin=Vout(:);
% % vout_test=kron(F*Da,kron(S,S))*vin
% % Vin=reshape(vin,Ms,N)
% V1=(Vout'.*da)';
% tt=[];
% for j=1:N
%     v=reshape(V1(:,j),M,M);
%     t1t=sqrt(2/(M+1))*dst(sqrt(2/(M+1))*dst(v).').';
%     tt(:,j)=t1t(:);
% end
% 
% Vout=1/sqrt(N)*(fft(tt.').');
% vout=Vout(:);
% 
% %% (5) (6) the fifth & the sixth step: compute vin=lam.^(-0.5).*vout & kron(inv(Da)*(F'),SS)*vin
% vin=lam.^(-0.5).*vout;
% V1=reshape(vin,Ms,N);
% tt=[];
% for j=1:N
%     v=reshape(V1(:,j),M,M);
%     t1t=sqrt(2/(M+1))*dst(sqrt(2/(M+1))*dst(v).').';
%     tt(:,j)=t1t(:);
% end
% Vout=sqrt(N)*(ifft(tt.').')*diag(1./da);
y=Vout(:);