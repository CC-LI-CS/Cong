%% laplacian with variable coefficients
clear;
clc;
x_l=0;
x_r=1;
y_l=0;
y_r=1;
dim=2;
levels=[6 8 10];
tol=1e-8;
maxit=1000;
DoF=[];
Iters=[];
t=[];
err=[];
for kkk=1:length(levels)
    level=levels(kkk);
    M=2^level-1; % space
    Ms=M*M;
    h=(x_r-x_l)/(M+1);
    fac=h^(-2);
    xgrid=(x_l:h:(x_r-h)).';
    xpgrid=xgrid+0.5*h;
    xgrid=xgrid(2:end);
    ygrid=(y_l:h:(y_r-h)).';
    ypgrid=ygrid+0.5*h;
    
    ygrid=ygrid(2:end);
    
    %% bar(beta)
    scalx=1/h^2; scaly=1/h^2;
    [Cxpmat,Cypmat,Cdmat]=getcoefmat(scalx,scaly,xgrid,xpgrid,ygrid,ypgrid);%get stepsize-scaled diffusion coefficients
    Cxpmat=-Cxpmat(2:end-1,:);
    Cypmat=-Cypmat(:,2:end-1);

    %% exact solution
    u_e=ones(M,M)+ones(M,M)*i;
    u0=u_e;
    Mx=M;
    My=M;
    %% get RHS
    %% rhs1=La*u_e
    rhs1=Cdmat.*u0+[Cxpmat;zeros(1,My)].*[u0(2:end,:);zeros(1,My)]+...
        [zeros(1,My);Cxpmat].*[zeros(1,My);u0(1:end-1,:)]+[Cypmat,zeros(Mx,1)].*[u0(:,2:end),zeros(Mx,1)]+...
        [zeros(Mx,1),Cypmat].*[zeros(Mx,1),u0(:,1:end-1)];

     

    %% Lambda = alpha1 + beta1*i
    alpha1=1;
    beta1=-100;
    lambda1=alpha1+beta1*i;
    %% rhs2=lambda1*u_e
    rhs2=lambda1*u_e;

    %% rhs
    rhs=rhs1+rhs2;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%test mgm

 y=level;

   
    %% fast
    real_ue=real(u_e);
    imag_ue=imag(u_e);
    real_rhs=real(rhs);
    imag_rhs=imag(rhs);

    %% rewitten linear system A1*z=f
    f=[imag_rhs(:);real_rhs(:)];
    DoF(kkk)=numel(f)
  
%% method 1 PRESB

tic;
[x,flag,relres,Iter] = gmres(@(vin)afun_var(alpha1,beta1,Cxpmat,Cypmat,Cdmat,2,M,Ms,1,vin),f,10,tol,maxit,@(tv)mgm_presb_tfun_var(Cxpmat,Cypmat,Cdmat,alpha1,beta1,tv,M,Ms,y,u_e));
x=x(1:Ms)+x(Ms+1:end)*i;
err(kkk)=norm(x-u_e(:),2)/norm(u_e(:),2)
Iter
%Iters(kkk)=Iter
t(kkk)=toc

end
