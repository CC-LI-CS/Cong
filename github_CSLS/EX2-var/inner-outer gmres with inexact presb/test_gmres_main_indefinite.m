%% laplacian with variable coefficients
clear;
clc;
x_l=0;
x_r=1;
y_l=0;
y_r=1;
dim=2;
levels=[6];
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
    Cxpmat=Cxpmat(2:end-1,:);
    Cypmat=Cypmat(:,2:end-1);

    %% exact solution
    u_e=randn(M,M)+randn(M,M)*1i;
    u0=u_e;
    Mx=M;
    My=M;
    %% get RHS
    %% rhs1=La*u_e
    rhs1=Cdmat.*u0+[-Cxpmat;zeros(1,My)].*[u0(2:end,:);zeros(1,My)]+...
        [zeros(1,My);-Cxpmat].*[zeros(1,My);u0(1:end-1,:)]+[-Cypmat,zeros(Mx,1)].*[u0(:,2:end),zeros(Mx,1)]+...
        [zeros(Mx,1),-Cypmat].*[zeros(Mx,1),u0(:,1:end-1)];

     

    %% Lambda = alpha1 + beta1*i
    alpha1=100;
    beta1=-100;
    lambda1=alpha1+beta1*1i;
    %% rhs2=lambda1*u_e
    rhs2=lambda1*u_e;

    %% rhs
    rhs=rhs1+rhs2;

    %% fast
    real_ue=real(u_e);
    imag_ue=imag(u_e);
    real_rhs=real(rhs);
    imag_rhs=imag(rhs);

    %% rewitten linear system A1*z=f
    f=[imag_rhs(:);real_rhs(:)];
    DoF(kkk)=numel(rhs(:))

% %% method 1 PRESB

tic;
offyseq=[-Cypmat,zeros(Mx,1)];
offyseq=offyseq(:);
offxseq=[-Cxpmat;zeros(1,My)];
offxseq=offxseq(:);
K=spdiags([offyseq,offxseq,Cdmat(:),[0;offxseq(1:end-1)],[zeros(Mx,1);offyseq(1:end-Mx)]],[-My,-1,0,1,My],Mx*My,Mx*My);
I=speye(size(K));
setup.type='nofill';  
eigB=beta1;

options.type = 'crout';
options.milu = 'row';
options.droptol = 0.1;
[L,U]=ilu( K+eigB*I ,options);

tic;
[x2,flag,relres,Iters2] = gmres(@(vin)afun_var_pmhss(alpha1,beta1,Cxpmat,Cypmat,Cdmat,2,M,Ms,1,vin),rhs(:),50,tol,maxit,@(tv)ilu_tfun_var(M,tv,Cxpmat,Cypmat,Cdmat,alpha1,beta1,Ms,L,U));
Iters_presb=Iters2
CPU_presb=toc
err2=norm(x2-u_e(:),2)/norm(u_e(:),2)
 
end
