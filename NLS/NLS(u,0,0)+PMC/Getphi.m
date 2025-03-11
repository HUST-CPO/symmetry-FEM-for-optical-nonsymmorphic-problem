function phi=Getphi(t1,t2,t3,u,v,w)
phi(1,1) = exp(1i*2*pi*(t1*u+t2*v+t3*w));
phi(1,2) = exp(1i*2*pi*(t1*u+t2*v+t3*w));