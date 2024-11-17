%% laplacian with variable coefficients
clear;
clc;
x_l=0;
x_r=1;
y_l=0;
y_r=1;
dim=2;
y=6:2:8;                                                             
slevels=[2.^y-1];
% levels=[6];
tol=1e-8;
maxit=1000;
DoF=[];
Iters=[];
t=[];
err=[];
for kkk=1:length(slevels)
    tic;  
    level=slevels(kkk);
    M=level; % space
    Ms=M*M;
    h=(x_r-x_l)/(M+1);
    fac=h^(-2);
    xgrid=(x_l:h:(x_r-h)).';
    xpgrid=xgrid+0.5*h;
    xgrid=xgrid(2:end);
    ygrid=(y_l:h:(y_r-h)).';
    ypgrid=ygrid+0.5*h;
    size(ypgrid)
    ygrid=ygrid(2:end);
    size(ygrid)
    %% bar(beta)
   
    scalx=1/h^2; scaly=1/h^2;
    [Cxpmat,Cypmat,Cdmat]=getcoefmat(scalx,scaly,xgrid,xpgrid,ygrid,ypgrid);%get stepsize-scaled diffusion coefficients
    Cxpmat=-Cxpmat(2:end-1,:);
    Cypmat=-Cypmat(:,2:end-1);

    barbeta=sqrt(min([Cxpmat(:);Cypmat(:)])*max([Cxpmat(:);Cypmat(:)]));
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

    size(rhs1)

    %% Lambda = alpha1 + beta1*i
    alpha1=1;
    beta1=-100;
    lambda1=alpha1+beta1*i;
    %% rhs2=lambda1*u_e
    rhs2=lambda1*u_e;

    %% rhs
    rhs=rhs1+rhs2;

    %% Finding kron(D,I)+kron(I,D) laplace
    theta=pi/(2*(M+1));
    seigs=sin((1:1:M)*theta).^2;
    sfacx=4;    sfacy=4;
    eigs=(kron(sfacx*ones(1,M),seigs)+kron(sfacy*seigs,ones(1,M)));
    eigs_P4=((barbeta*eigs'+alpha1.*ones(Ms,1)).^2+beta1^2.*ones(Ms,1)).^(1/2);
    eigs_P4=[eigs_P4;eigs_P4];
    % cT=tcoeff
    %% T
    %% fast
    real_ue=real(u_e);
    imag_ue=imag(u_e);
    real_rhs=real(rhs);
    imag_rhs=imag(rhs);

    %% rewitten linear system A1*z=f
    f=[imag_rhs(:);real_rhs(:)];
    DoF(kkk)=numel(f);

    %%
    [x,flag,relres,Iters(kkk)] = minres(@(vin)afun_var(alpha1,beta1,Cxpmat,Cypmat,Cdmat,2,M,Ms,1,vin) ,f,tol,maxit,@(tv)tfun_var(eigs_P4,tv,1,M,Ms));
    x=x(1:Ms)+x(Ms+1:end)*i;
    t(kkk)=toc;
    size(x)
    size(u_e)
    err(kkk)=norm(x-u_e(:),2)/norm(u_e(:),2)

end
DoF
Iters
t
err;


