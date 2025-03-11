function phi=Getglidephi(u,tage)
if tage==1
    phi(1,1) = exp(1i*pi*u);
    phi(1,2) = exp(1i*pi*(1+u));
elseif tage==2
    phi(1,1) = 1;
    phi(1,2) = -1;
elseif tage==3
    phi(1,1) = 1;
    phi(1,2) = exp(1i*2*pi*u);
elseif tage==4
    phi(1,1) = exp(1i*pi*u);
    phi(1,2) = exp(1i*pi*(1+u));
elseif tage==5
    phi(1,1) = exp(1i*pi*u);
    phi(1,2) = exp(1i*pi*(1+u));
elseif tage==6
    phi(1,1) = 1;
    phi(1,2) = 1;
end