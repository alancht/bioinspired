function [G]=three_stokeslets_3D(x,y,z,xi,yi,zi,mu)

%mu=1;

dx=x-xi;
dy=y-yi;
dz=z-zi;
r=sqrt(dx^2+dy^2+dz^2);

Gxx=1/8/pi/mu*(1/r+dx*dx/r^3);
Gyy=1/8/pi/mu*(1/r+dy*dy/r^3);
Gzz=1/8/pi/mu*(1/r+dz*dz/r^3);

Gxy=1/8/pi/mu*(0/r+(dx*dy)/r^3);
Gxz=1/8/pi/mu*(0/r+(dx*dz)/r^3);
Gyz=1/8/pi/mu*(0/r+(dy*dz)/r^3);

Gyx=Gxy;
Gzx=Gxz;
Gzy=Gyz;

G=[Gxx,Gxy,Gxz;Gyx,Gyy,Gyz;Gzx,Gzy,Gzz];

end