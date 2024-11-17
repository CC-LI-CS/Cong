function y = ilu_tfun_var(M,tv,Cxpmat,Cypmat,Cdmat,alpha1,beta1,Ms,L,U)
real_tv=real(tv);
imag_tv=imag(tv);
z1=[imag_tv(:);real_tv(:)];
x = gmres(@(vin)afun_var(0,beta1,Cxpmat,Cypmat,Cdmat,2,M,Ms,1,vin),z1,60,1e-9,10000,@(ttv)ilu_presb_tfun_var(M,ttv,Cxpmat,Cypmat,Cdmat,0,Ms,L,U));
x=x(1:Ms)+x(Ms+1:end)*i;
y=(beta1/(abs(alpha1)*1i+beta1))*x;
