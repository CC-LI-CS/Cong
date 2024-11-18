function y = mgm_presb_tfun_var(Cxpmat,Cypmat,Cdmat,alpha1,beta1,tv,M,Ms,level,u_e)

Tv=tv(1:Ms);
Tv1=tv(Ms+1:end);
y2_1=Tv1-Tv;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha1=-100;beta1=25;
eigB=alpha1+beta1;
z1=mgm_inv(y2_1(:),eigB,level,u_e(:));
z1=-z1;
%% rhs1=La*u_e
u0=reshape(z1,M,M);
My=M;Mx=M;
z11=Cdmat.*u0+[Cxpmat;zeros(1,My)].*[u0(2:end,:);zeros(1,My)]+...
        [zeros(1,My);Cxpmat].*[zeros(1,My);u0(1:end-1,:)]+[Cypmat,zeros(Mx,1)].*[u0(:,2:end),zeros(Mx,1)]+...
        [zeros(Mx,1),Cypmat].*[zeros(Mx,1),u0(:,1:end-1)];

z2=Tv-z11(:)-alpha1*z1;
z3=mgm_inv(z2(:),alpha1,level,u_e(:));

y=[z3;z1-z3];
