function [Cxpmat,Cypmat,Cdmat]=getcoefmat(scalx,scaly,xgrid,xpgrid,ygrid,ypgrid)
xgrid=xgrid(:);
xpgrid=xpgrid(:);
ygrid=ygrid(:);
ypgrid=ypgrid(:);
% length(ygrid)
% length(xpgrid)
[Ycoords,Xpcoords]=meshgrid(ygrid,xpgrid);
Cxpmat=scalx*diffu_coeff(Xpcoords,Ycoords);

[Ypcoords,Xcoords]=meshgrid(ypgrid,xgrid);
Cypmat=scaly*diffu_coeff(Xcoords,Ypcoords);

% Cxpmat=fac*(40*ones(size((kron(ygrid',xpgrid))))+kron(ones(size(ygrid')),xpgrid.^3.5)+kron(ygrid'.^3.5,ones(size(xpgrid))));
% Cypmat=fac*(40*ones(size((kron(ypgrid',xgrid))))+kron(ones(size(ypgrid')),xgrid.^3.5)+kron(ypgrid'.^3.5,ones(size(xgrid))));
% Cdmat=Cxpmat(1:end-1,:)+Cxpmat(2:end,:)+Cypmat(:,1:end-1)+Cypmat(:,2:end);

% size(Cypmat)
% size(Cxpmat)
Cdmat=Cxpmat(1:end-1,:)+Cxpmat(2:end,:)+Cypmat(:,1:end-1)+Cypmat(:,2:end);
end