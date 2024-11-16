function vout=dstn(vin,dim)
%discrete sine transform along a dimension of vin specified by dim
%vin must be real
siza=size(vin);
fillsiz1=siza;
fillsiz1(dim)=1;
fillsiz2=siza;
fillsiz2(dim)=siza(dim)+1;
nd=length(siza);
truccomline='vout=vout(';
for ii=1:nd
    if ii~=dim
       truccomline=strcat(truccomline,':'); 
    else
       truccomline=strcat(truccomline,'2:',num2str(siza(dim)),'+1'); 
    end
    if ii~=nd
        truccomline=strcat(truccomline,',');
    end
end
truccomline=strcat(truccomline,');');
vout=imag(fft(cat(dim,zeros(fillsiz1),-vin,zeros(fillsiz2)),[],dim));
eval(truccomline)
end