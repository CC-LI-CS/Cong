function y = ilu_presb_tfun_var(M,ttv,Cxpmat,Cypmat,Cdmat,alpha1,Ms,L,U)

y1=ttv(1:Ms);
y2=ttv(Ms+1:end)-y1;
Z1=-U\(L\y2(:));
%% rhs1=La*u_e
temp=reshape(Z1,M,M);
My=M;Mx=M;
temp=Cdmat.*temp+[-Cxpmat;zeros(1,My)].*[temp(2:end,:);zeros(1,My)]+...
        [zeros(1,My);-Cxpmat].*[zeros(1,My);temp(1:end-1,:)]+[-Cypmat,zeros(Mx,1)].*[temp(:,2:end),zeros(Mx,1)]+...
        [zeros(Mx,1),-Cypmat].*[zeros(Mx,1),temp(:,1:end-1)];

temp=temp(:)+alpha1*Z1(:);
Z2=y1-temp(:);

% z3=mgm_inv(z2(:),alpha1,level,u_e(:));
Z3=U\(L\Z2(:));

y=[Z3;Z1-Z3];
