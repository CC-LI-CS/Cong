%% laplacian with constant coefficients
clear;
clc;
x_l=0;
x_r=1;
y_l=0;
y_r=1;
N=1;
slevels=2^7-1;
tol=1e-8;
maxit=1000;
DoF=[];
Iters=[];
err=[];
t=[];
for kkk=1:length(slevels)
    level=slevels(kkk);
    M=level; % space
    Ms=M*M
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
    %% Finding kron(D,I)+kron(I,D) laplace
    theta=pi/(2*(M+1));
    seigs=sin((1:1:M)*theta).^2;
    sfacx=4;    sfacy=4;
    eigs=(kron(sfacx*ones(1,M),seigs)+kron(sfacy*seigs,ones(1,M)));

    %% L1=dt*hx^(-2)*(kron(W,IM)+kron(IM,W));
    alpha1=100;
    beta1=-100;
    lambda1=alpha1+beta1*i;
    eigs_P4=((h^(-2)*eigs'+alpha1.*ones(Ms,1)).^2+beta1^2.*ones(Ms,1)).^(1/2);
    eigs_P4=[eigs_P4;eigs_P4];
    %% rhs
    u_exact=ones(Ms,1)+ones(Ms,1)*i;
    Nt=1;
    dt=1;
    Mx=M;
    My=M;
    scalx=dt/h^2; scaly=dt/h^2;
    vin=u_exact(:);
    vin=reshape(vin,Nt,Mx,My);
    spatialpart=scalx*[2*vin(:,1,:)-vin(:,2,:),2*vin(:,2:end-1,:)-vin(:,3:end,:)-vin(:,1:end-2,:),2*vin(:,end,:)-vin(:,end-1,:)]+...
        scaly*cat(3,2*vin(:,:,1)-vin(:,:,2),2*vin(:,:,2:end-1)-vin(:,:,3:end)-vin(:,:,1:end-2),2*vin(:,:,end)-vin(:,:,end-1));

    rhs1=spatialpart(:);
    %% rhs1=L*u_e
    %% rhs2=lambda1*u_e
    rhs2=lambda1*u_exact;
    %% rhs
    rhs=rhs1+rhs2;

    %% fast
    real_ue=real(u_exact);
    imag_ue=imag(u_exact);
    real_rhs=real(rhs);
    imag_rhs=imag(rhs);
    f=[imag_rhs(:);real_rhs(:)];
    DoF(kkk)=numel(f)/2;
    tic;
    [x,flag,relres,Iters(kkk)] = minres(@(vin)afun_con(alpha1,beta1,scalx,scaly,1,Mx,My,vin),f,tol,maxit,@(tv)tfun_con(eigs_P4,tv,1,M,Ms));
    x=x(1:Ms)+x(Ms+1:end)*i;
    t(kkk)=toc;
    size(x)
    err(kkk)=norm(x-u_exact,2)/norm(u_exact,2)
end
DoF
Iters
t
err
