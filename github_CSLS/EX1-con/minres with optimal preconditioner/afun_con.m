function vout = afun_con(alpha1,beta1,scalx,scaly,Nt,Mx,My,vin)
% this is for computing YT times a vector vin and
%returning the result in vout;
%Nt is the number of temporal grid points
%Ns is the number of spatial grid points 
%scalx=1/(hx^2),scaly=1/(hy^2), ax and ay are the constant
%diffusion coefficients, Mx (My, respectively) is the number of spatial grid points
% along x (y, respectively) direction.

%computing -\nabda times vin
Ms=Mx*My;
% vin1=reshape(vin(1:Ms),Nt,Mx,My);
% spatialpart=scalx*[2*vin(:,1,:)-vin(:,2,:),2*vin(:,2:end-1,:)-vin(:,3:end,:)-vin(:,1:end-2,:),2*vin(:,end,:)-vin(:,end-1,:)]+...
%             scaly*cat(3,2*vin(:,:,1)-vin(:,:,2),2*vin(:,:,2:end-1)-vin(:,:,3:end)-vin(:,:,1:end-2),2*vin(:,:,end)-vin(:,:,end-1));
% 
% vout=spatialpart(:);
            
%% 
x1=vin(1:Ms);
x2=vin(Ms+1:end);

%% L*x1
vin1=reshape(x1,Nt,Mx,My);
vout1=scalx*[2*vin1(:,1,:)-vin1(:,2,:),2*vin1(:,2:end-1,:)-vin1(:,3:end,:)-vin1(:,1:end-2,:),2*vin1(:,end,:)-vin1(:,end-1,:)]+...
            scaly*cat(3,2*vin1(:,:,1)-vin1(:,:,2),2*vin1(:,:,2:end-1)-vin1(:,:,3:end)-vin1(:,:,1:end-2),2*vin1(:,:,end)-vin1(:,:,end-1));
vout1=vout1(:);

%% L*x2
vin2=reshape(x2,Nt,Mx,My);                                    
vout2=scalx*[2*vin2(:,1,:)-vin2(:,2,:),2*vin2(:,2:end-1,:)-vin2(:,3:end,:)-vin2(:,1:end-2,:),2*vin2(:,end,:)-vin2(:,end-1,:)]+...
            scaly*cat(3,2*vin2(:,:,1)-vin2(:,:,2),2*vin2(:,:,2:end-1)-vin2(:,:,3:end)-vin2(:,:,1:end-2),2*vin2(:,:,end)-vin2(:,:,end-1));
vout2=vout2(:);

%% 
vout=[beta1*x1+alpha1*x2+vout2;alpha1*x1-beta1*x2+vout1];
end