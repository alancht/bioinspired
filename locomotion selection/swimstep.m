function [r0_l,phi1,phi2,dr0_l,deye_l,dyaw_l,dp_l,ex_b,ey_b,ez_b,R1ns,R2ns,Work,w1all,w2all,W_motor,info,d_ff,V0]=swimstep(r0_l,phi1,phi2,dr0_l,deye_l,dyaw_l,dp_l,ex_b,ey_b,ez_b,R1ns,R2ns,l1p,l2p,h1p,h2p,nt,dt)

ex_l=[1;0;0];          % coordinate axis in lab view
ey_l=[0;1;0];
ez_l=[0;0;1];

theta1=180/180*pi;
theta2=210/180*pi;

beta1=0.3;
beta2=0.3;

c10=0.7;
c20=0.7;
c1=c10;
c2=c20;
%Flagella waveforms
b1 = 0.8;
b2 = 1.5;

l10=l1p*b1;   %typical cell size    trans
l20=l2p*b1;   %cis
h10=h1p*b1;
h20=h2p*b1;
l1=l10;
l2=l20;
h1=h10;
h2=h20;

% light direction and strength
Ix=0;
Iy=0;
Iz=0;

I0=0;
p1=-0.1;       % 
p2=+0.1; 

Ra10=0.1572*b2;
Rb10=0.2803*b2;
Ra20=0.1649*b2;
Rb20=0.2906*b2;


Ra1=Ra10;
Rb1=Rb10;
Ra2=Ra20;
Rb2=Rb20;

a0=0.5;    %radius of body sphere
a1=0.075;     %radius of flagella sphere
a2=0.075;    
mu=(10.6/10)^2;

%If ellipse
phi11=atan2(sin(phi1)*Ra1,cos(phi1)*Rb1);
phi22=atan2(sin(phi2)*Ra2,cos(phi2)*Rb2);
R1=sqrt(Ra1*Ra1*cos(phi11)*cos(phi11) + Rb1*Rb1*sin(phi11)*sin(phi11));
R2=sqrt(Ra2*Ra2*cos(phi22)*cos(phi22) + Rb2*Rb2*sin(phi22)*sin(phi22));
  
%Recover the useful parameters
a_orbit1 = pi/2-beta1;
v_orbit1 = rotate(ez_b,ex_b,a_orbit1);
a_orbit2 = beta2;
v_orbit2 = rotate(ey_b,ex_b,a_orbit2);
vn1_0 = ex_b;
vn2_0 = -ex_b;

F1n_l = rotate(vn1_0,v_orbit1,phi1);
F2n_l = rotate(vn2_0,v_orbit2,phi2);
F1t_l = rotate(F1n_l,v_orbit1,pi/2); %tangent direction
F2t_l = rotate(F2n_l,v_orbit2,pi/2);

r01 = -l1*ex_b + h1*ez_b;
r02 = l2*ex_b + h2*ez_b;

r1_l = r0_l + r01 + R1*F1n_l;
r2_l = r0_l + r02 + R2*F2n_l;

r1_lp = r0_l+r01/b1+R1/b2*F1n_l;
r2_lp = r0_l+r02/b1+R2/b2*F2n_l;

F1t=1+c10*cos(phi1+theta1);
F2t=1+c20*cos(phi2+theta2);

%% define data matrix
aphi1n=zeros(1,nt);
aphi2n=zeros(1,nt);
adphi1n=zeros(1,nt);
adphi2n=zeros(1,nt);
awn=zeros(3,nt);
ar0n=zeros(3,nt);
ar1n=zeros(3,nt);
ar2n=zeros(3,nt);
ar1np=zeros(3,nt);
ar2np=zeros(3,nt);
adr0n=zeros(3,nt);
aden=zeros(3,nt);
adyn=zeros(3,nt);
adpn=zeros(3,nt);
aorbit1=zeros(3,nt);
aorbit2=zeros(3,nt);
avn1_0=zeros(3,nt);
avn2_0=zeros(3,nt);
ar01=zeros(3,nt);
ar02=zeros(3,nt);
aF1n_l=zeros(3,nt);
aF1t_l=zeros(3,nt);
aF2n_l=zeros(3,nt);
aF2t_l=zeros(3,nt);
aI1=zeros(1,nt);
apI1=zeros(1,nt);
ac1=zeros(1,nt);
ac2=zeros(1,nt);
ah1=zeros(1,nt);
ah2=zeros(1,nt);
aR1=zeros(1,nt);
aR2=zeros(1,nt);
abwn=zeros(3,nt);
abT1n=zeros(3,nt);
abT2n=zeros(3,nt);
aF1n=zeros(1,nt);
aF2n=zeros(1,nt);
wfla = [];
Work = 0;
w1all = 0;
w2all = 0;
W_motor = 0;
W_fla = 0;
d_ff=[];
wleft=[];
wright=[];
V0=[];
ang1 = [];
ang2 = [];
for i=1:nt

    [X,A,B]=Oseen_tensor(r0_l,r1_l,r2_l,phi1,phi2,theta1,theta2,F1n_l,F1t_l,F2n_l,F2t_l,v_orbit1,v_orbit2,mu,Ra1,Rb1,Ra2,Rb2,R1ns,R2ns,a1,a2,a0,c1,c2);

    F1n=X(1);
    F2n=X(2);
    wx=X(3);
    wy=X(4);
    wz=X(5);
    dphi1=X(6);
    dphi2=X(7);
    F0x=X(8);
    F0y=X(9);
    F0z=X(10);
    v1n=X(11);
    v1t=X(12);
    v2n=X(13);
    v2t=X(14);
    v0x=X(15);
    v0y=X(16);
    v0z=X(17);
    F1p=X(18);
    F2p=X(19);
    v1p=X(20);
    v2p=X(21);
    F1t=1+c1*cos(phi1+theta1);
    F2t=1+c2*cos(phi2+theta2);

    r0_old = r0_l;

%% calculate torque
   w=[wx; wy; wz];
   r10=r1_l-r0_l;
   r20=r2_l-r0_l;
   F0=[F0x; F0y; F0z];
   F1=F1n*F1n_l + F1t*F1t_l + F1p*v_orbit1;
   F2=F2n*F2n_l + F2t*F2t_l + F2p*v_orbit2;

   vleft7 = 8*pi*mu*a0*a0*a0*w + cross(r10,F1) + cross(r20,F2);
   %torq_old = cross(r0_old,F_old) + 8*pi*mu*a0*a0*a0*w;

    T1 = cross(r10,F1);
    T2 = cross(r20,F2);
    


    Ldr = dr0_l;
    Ldy = dyaw_l;
    Ldp = dp_l;
    L_M = [Ldr(1), Ldy(1), Ldp(1);  Ldr(2), Ldy(2), Ldp(2);  Ldr(3), Ldy(3), Ldp(3)];

    Bw = L_M \ w;
    BT1 = L_M \ T1;
    BT2 = L_M \ T2;

    abwn(:,i)=Bw;
    abT1n(:,i)=BT1;
    abT2n(:,i)=BT2;



%% resonding to light

    I1=-(Ix*deye_l(1)+Iy*deye_l(2)+Iz*deye_l(3));    % evec: direction of maximal light sensitivity
        
    if (Ix==0 && Iy==0 && Iz==0)
            p=0;      
    end

    I1=1+I0*I1*heaviside(I1);

    h1=h10+p1*log(I1);
    h2=h20+p2*log(I1);
% 
%     l1=l10-p1*log(I1);
%     l2=l20-p2*log(I1);

    Ra1=Ra10+0.4*p1*log(I1);
    Ra2=Ra20+0.4*p2*log(I1);
% 
    Rb1=Rb10+0.8*p1*log(I1);
    Rb2=Rb20+0.8*p2*log(I1);

%% Update position
    aphi1n(:,i)=phi1;
    aphi2n(:,i)=phi2;
    adphi1n(:,i)=dphi1;
    adphi2n(:,i)=dphi2;
    awn(:,i)=[X(3);X(4);X(5)];
    ar0n(:,i)=r0_l;     
    ar1n(:,i)=r1_l;
    ar2n(:,i)=r2_l;
    ar1np(:,i)=r1_lp;
    ar2np(:,i)=r2_lp;
    adr0n(:,i)=dr0_l;
    aden(:,i)=deye_l;
    adyn(:,i)=dyaw_l;
    adpn(:,i)=dp_l;
    aorbit1(:,i)=v_orbit1;
    aorbit2(:,i)=v_orbit2;
    avn1_0(:,i)=vn1_0;
    avn2_0(:,i)=vn2_0;
    ar01(:,i)=r01;
    ar02(:,i)=r02;
    aF1n_l(:,i)=F1n_l;
    aF1t_l(:,i)=F1t_l;
    aF2n_l(:,i)=F2n_l;
    aF2t_l(:,i)=F2t_l;
    aI1(:,i)=I1;
    apI1(:,i)=log(I1);
    ac1(:,i)=c1;
    ac2(:,i)=c2;
    ah1(:,i)=h1;
    ah2(:,i)=h2;
    aR1(:,i)=(Ra1+Rb1)/2;
    aR2(:,i)=(Ra2+Rb2)/2;
    aF1n(:,i)=F1n;
    aF2n(:,i)=F2n;    

%% RK4
% k1
    v0x_k1 = X(15);
    v0y_k1 = X(16);
    v0z_k1 = X(17);
    v1x_k1 = X(11)*F1n_l(1) + X(12)*F1t_l(1) + X(20)*v_orbit1(1);
    v1y_k1 = X(11)*F1n_l(2) + X(12)*F1t_l(2) + X(20)*v_orbit1(2);
    v1z_k1 = X(11)*F1n_l(3) + X(12)*F1t_l(3) + X(20)*v_orbit1(3);
    v2x_k1 = X(13)*F2n_l(1) + X(14)*F2t_l(1) + X(21)*v_orbit2(1);
    v2y_k1 = X(13)*F2n_l(2) + X(14)*F2t_l(2) + X(21)*v_orbit2(2);
    v2z_k1 = X(13)*F2n_l(3) + X(14)*F2t_l(3) + X(21)*v_orbit2(3);

    F1_k1=X(1)*F1n_l + F1t*F1t_l + X(19)*v_orbit1;
    F2_k1=X(2)*F2n_l + F2t*F2t_l + X(20)*v_orbit2;

    F0_k1 = [X(8),X(9),X(10)]';

    dphi1_k1 = X(6);
    dphi2_k1 = X(7);
    w0x_k1 = X(3);
    w0y_k1 = X(4);
    w0z_k1 = X(5);

 % k2
    r0_l_k2 = r0_l;
    r1_l_k2 = r1_l;
    r2_l_k2 = r2_l;
    phi1_k2 = phi1;
    phi2_k2 = phi2;
    v_orbit1_k2 = v_orbit1;
    v_orbit2_k2 = v_orbit2;
    vn1_0_k2 = vn1_0;
    vn2_0_k2 = vn2_0;


    r0_l_k2(1)=r0_l_k2(1)+dt/2*v0x_k1;
    r0_l_k2(2)=r0_l_k2(2)+dt/2*v0y_k1;
    r0_l_k2(3)=r0_l_k2(3)+dt/2*v0z_k1;

    r1_l_k2(1)=r1_l_k2(1)+dt/2*v1x_k1;
    r1_l_k2(2)=r1_l_k2(2)+dt/2*v1y_k1;
    r1_l_k2(3)=r1_l_k2(3)+dt/2*v1z_k1;

    r2_l_k2(1)=r2_l_k2(1)+dt/2*v2x_k1;
    r2_l_k2(2)=r2_l_k2(2)+dt/2*v2y_k1;
    r2_l_k2(3)=r2_l_k2(3)+dt/2*v2z_k1;

    phi1_k2=phi1_k2+dt/2*dphi1_k1;
    phi2_k2=phi2_k2+dt/2*dphi2_k1;

    v_orbit1_k2     =  rotate(v_orbit1_k2,ex_l,w0x_k1*dt/2); 
    v_orbit2_k2     =  rotate(v_orbit2_k2,ex_l,w0x_k1*dt/2); 
    vn1_0_k2        =  rotate(vn1_0_k2,ex_l,w0x_k1*dt/2);
    vn2_0_k2        =  rotate(vn2_0_k2,ex_l,w0x_k1*dt/2);

    v_orbit1_k2     =  rotate(v_orbit1_k2,ey_l,w0y_k1*dt/2); 
    v_orbit2_k2     =  rotate(v_orbit2_k2,ey_l,w0y_k1*dt/2); 
    vn1_0_k2        =  rotate(vn1_0_k2,ey_l,w0y_k1*dt/2);
    vn2_0_k2        =  rotate(vn2_0_k2,ey_l,w0y_k1*dt/2);

    v_orbit1_k2     =  rotate(v_orbit1_k2,ez_l,w0z_k1*dt/2); 
    v_orbit2_k2     =  rotate(v_orbit2_k2,ez_l,w0z_k1*dt/2); 
    vn1_0_k2        =  rotate(vn1_0_k2,ez_l,w0z_k1*dt/2);
    vn2_0_k2        =  rotate(vn2_0_k2,ez_l,w0z_k1*dt/2);

    F1n_l_k2        =  rotate(vn1_0_k2,v_orbit1_k2,phi1_k2);
    F2n_l_k2        =  rotate(vn2_0_k2,v_orbit2_k2,phi2_k2);
    F1t_l_k2        =  rotate(F1n_l_k2,v_orbit1_k2,pi/2);
    F2t_l_k2        =  rotate(F2n_l_k2,v_orbit2_k2,pi/2);

    [X_k2]=Oseen_tensor(r0_l_k2,r1_l_k2,r2_l_k2,phi1_k2,phi2_k2,theta1,theta2,F1n_l_k2,F1t_l_k2,F2n_l_k2,F2t_l_k2,v_orbit1_k2,v_orbit2_k2,mu,Ra1,Rb1,Ra2,Rb2,R1ns,R2ns,a1,a2,a0,c1,c2);

    v0x_k2 = X_k2(15);
    v0y_k2 = X_k2(16);
    v0z_k2 = X_k2(17);
    v1x_k2 = X_k2(11)*F1n_l_k2(1) + X_k2(12)*F1t_l_k2(1) + X_k2(20)*v_orbit1_k2(1);
    v1y_k2 = X_k2(11)*F1n_l_k2(2) + X_k2(12)*F1t_l_k2(2) + X_k2(20)*v_orbit1_k2(2);
    v1z_k2 = X_k2(11)*F1n_l_k2(3) + X_k2(12)*F1t_l_k2(3) + X_k2(20)*v_orbit1_k2(3);
    v2x_k2 = X_k2(13)*F2n_l_k2(1) + X_k2(14)*F2t_l_k2(1) + X_k2(21)*v_orbit2_k2(1);
    v2y_k2 = X_k2(13)*F2n_l_k2(2) + X_k2(14)*F2t_l_k2(2) + X_k2(21)*v_orbit2_k2(2);
    v2z_k2 = X_k2(13)*F2n_l_k2(3) + X_k2(14)*F2t_l_k2(3) + X_k2(21)*v_orbit2_k2(3);

    F1_k2=X_k2(1)*F1n_l_k2 + F1t*F1t_l_k2 + X_k2(19)*v_orbit1_k2;
    F2_k2=X_k2(2)*F2n_l_k2 + F2t*F2t_l_k2 + X_k2(20)*v_orbit2_k2;
    F0_k2 = [X_k2(8),X_k2(9),X_k2(10)]';

    dphi1_k2 = X_k2(6);
    dphi2_k2 = X_k2(7);
    w0x_k2 = X_k2(3);
    w0y_k2 = X_k2(4);
    w0z_k2 = X_k2(5);
 
 % k3
    r0_l_k3 = r0_l;
    r1_l_k3 = r1_l;
    r2_l_k3 = r2_l;
    phi1_k3 = phi1;
    phi2_k3 = phi2;
    v_orbit1_k3 = v_orbit1;
    v_orbit2_k3 = v_orbit2;
    vn1_0_k3 = vn1_0;
    vn2_0_k3 = vn2_0;


    r0_l_k3(1)=r0_l_k3(1)+dt/2*v0x_k2;
    r0_l_k3(2)=r0_l_k3(2)+dt/2*v0y_k2;
    r0_l_k3(3)=r0_l_k3(3)+dt/2*v0z_k2;

    r1_l_k3(1)=r1_l_k3(1)+dt/2*v1x_k2;
    r1_l_k3(2)=r1_l_k3(2)+dt/2*v1y_k2;
    r1_l_k3(3)=r1_l_k3(3)+dt/2*v1z_k2;

    r2_l_k3(1)=r2_l_k3(1)+dt/2*v2x_k2;
    r2_l_k3(2)=r2_l_k3(2)+dt/2*v2y_k2;
    r2_l_k3(3)=r2_l_k3(3)+dt/2*v2z_k2;

    phi1_k3=phi1_k3+dt/2*dphi1_k2;
    phi2_k3=phi2_k3+dt/2*dphi2_k2;

    v_orbit1_k3     =  rotate(v_orbit1_k3,ex_l,w0x_k2*dt/2); 
    v_orbit2_k3     =  rotate(v_orbit2_k3,ex_l,w0x_k2*dt/2); 
    vn1_0_k3        =  rotate(vn1_0_k3,ex_l,w0x_k2*dt/2);
    vn2_0_k3        =  rotate(vn2_0_k3,ex_l,w0x_k2*dt/2);

    v_orbit1_k3     =  rotate(v_orbit1_k3,ey_l,w0y_k2*dt/2); 
    v_orbit2_k3     =  rotate(v_orbit2_k3,ey_l,w0y_k2*dt/2); 
    vn1_0_k3        =  rotate(vn1_0_k3,ey_l,w0y_k2*dt/2);
    vn2_0_k3        =  rotate(vn2_0_k3,ey_l,w0y_k2*dt/2);

    v_orbit1_k3     =  rotate(v_orbit1_k3,ez_l,w0z_k2*dt/2); 
    v_orbit2_k3     =  rotate(v_orbit2_k3,ez_l,w0z_k2*dt/2); 
    vn1_0_k3        =  rotate(vn1_0_k3,ez_l,w0z_k2*dt/2);
    vn2_0_k3        =  rotate(vn2_0_k3,ez_l,w0z_k2*dt/2);

    F1n_l_k3        =  rotate(vn1_0_k3,v_orbit1_k3,phi1_k3);
    F2n_l_k3        =  rotate(vn2_0_k3,v_orbit2_k3,phi2_k3);
    F1t_l_k3        =  rotate(F1n_l_k3,v_orbit1_k3,pi/2);
    F2t_l_k3        =  rotate(F2n_l_k3,v_orbit2_k3,pi/2);

    [X_k3]=Oseen_tensor(r0_l_k3,r1_l_k3,r2_l_k3,phi1_k3,phi2_k3,theta1,theta2,F1n_l_k3,F1t_l_k3,F2n_l_k3,F2t_l_k3,v_orbit1_k3,v_orbit2_k3,mu,Ra1,Rb1,Ra2,Rb2,R1ns,R2ns,a1,a2,a0,c1,c2);

    v0x_k3 = X_k3(15);
    v0y_k3 = X_k3(16);
    v0z_k3 = X_k3(17);
    v1x_k3 = X_k3(11)*F1n_l_k3(1) + X_k3(12)*F1t_l_k3(1) + X_k3(20)*v_orbit1_k3(1);
    v1y_k3 = X_k3(11)*F1n_l_k3(2) + X_k3(12)*F1t_l_k3(2) + X_k3(20)*v_orbit1_k3(2);
    v1z_k3 = X_k3(11)*F1n_l_k3(3) + X_k3(12)*F1t_l_k3(3) + X_k3(20)*v_orbit1_k3(3);
    v2x_k3 = X_k3(13)*F2n_l_k3(1) + X_k3(14)*F2t_l_k3(1) + X_k3(21)*v_orbit2_k3(1);
    v2y_k3 = X_k3(13)*F2n_l_k3(2) + X_k3(14)*F2t_l_k3(2) + X_k3(21)*v_orbit2_k3(2);
    v2z_k3 = X_k3(13)*F2n_l_k3(3) + X_k3(14)*F2t_l_k3(3) + X_k3(21)*v_orbit2_k3(3);
    
    F1_k3=X_k3(1)*F1n_l_k3 + F1t*F1t_l_k3 + X_k3(19)*v_orbit1_k3;
    F2_k3=X_k3(2)*F2n_l_k3 + F2t*F2t_l_k3 + X_k3(20)*v_orbit2_k3;
    F0_k3 = [X_k3(8),X_k3(9),X_k3(10)]';
    dphi1_k3 = X_k3(6);
    dphi2_k3 = X_k3(7);
    w0x_k3 = X_k3(3);
    w0y_k3 = X_k3(4);
    w0z_k3 = X_k3(5);

 % k4
    r0_l_k4 = r0_l;
    r1_l_k4 = r1_l;
    r2_l_k4 = r2_l;
    phi1_k4 = phi1;
    phi2_k4 = phi2;
    v_orbit1_k4 = v_orbit1;
    v_orbit2_k4 = v_orbit2;
    vn1_0_k4 = vn1_0;
    vn2_0_k4 = vn2_0;


    r0_l_k4(1)=r0_l_k4(1)+dt*v0x_k3;
    r0_l_k4(2)=r0_l_k4(2)+dt*v0y_k3;
    r0_l_k4(3)=r0_l_k4(3)+dt*v0z_k3;

    r1_l_k4(1)=r1_l_k4(1)+dt*v1x_k3;
    r1_l_k4(2)=r1_l_k4(2)+dt*v1y_k3;
    r1_l_k4(3)=r1_l_k4(3)+dt*v1z_k3;

    r2_l_k4(1)=r2_l_k4(1)+dt*v2x_k3;
    r2_l_k4(2)=r2_l_k4(2)+dt*v2y_k3;
    r2_l_k4(3)=r2_l_k4(3)+dt*v2z_k3;

    phi1_k4=phi1_k4+dt*dphi1_k3;
    phi2_k4=phi2_k4+dt*dphi2_k3;

    v_orbit1_k4     =  rotate(v_orbit1_k4,ex_l,w0x_k3*dt); 
    v_orbit2_k4     =  rotate(v_orbit2_k4,ex_l,w0x_k3*dt); 
    vn1_0_k4        =  rotate(vn1_0_k4,ex_l,w0x_k3*dt);
    vn2_0_k4        =  rotate(vn2_0_k4,ex_l,w0x_k3*dt);

    v_orbit1_k4     =  rotate(v_orbit1_k4,ey_l,w0y_k3*dt); 
    v_orbit2_k4     =  rotate(v_orbit2_k4,ey_l,w0y_k3*dt); 
    vn1_0_k4        =  rotate(vn1_0_k4,ey_l,w0y_k3*dt);
    vn2_0_k4        =  rotate(vn2_0_k4,ey_l,w0y_k3*dt);

    v_orbit1_k4     =  rotate(v_orbit1_k4,ez_l,w0z_k3*dt); 
    v_orbit2_k4     =  rotate(v_orbit2_k4,ez_l,w0z_k3*dt); 
    vn1_0_k4        =  rotate(vn1_0_k4,ez_l,w0z_k3*dt);
    vn2_0_k4        =  rotate(vn2_0_k4,ez_l,w0z_k3*dt);

    F1n_l_k4        =  rotate(vn1_0_k4,v_orbit1_k4,phi1_k4);
    F2n_l_k4        =  rotate(vn2_0_k4,v_orbit2_k4,phi2_k4);
    F1t_l_k4        =  rotate(F1n_l_k4,v_orbit1_k4,pi/2);
    F2t_l_k4        =  rotate(F2n_l_k4,v_orbit2_k4,pi/2);

    [X_k4]=Oseen_tensor(r0_l_k4,r1_l_k4,r2_l_k4,phi1_k4,phi2_k4,theta1,theta2,F1n_l_k4,F1t_l_k4,F2n_l_k4,F2t_l_k4,v_orbit1_k4,v_orbit2_k4,mu,Ra1,Rb1,Ra2,Rb2,R1ns,R2ns,a1,a2,a0,c1,c2);

    v0x_k4 = X_k4(15);
    v0y_k4 = X_k4(16);
    v0z_k4 = X_k4(17);
    v1x_k4 = X_k4(11)*F1n_l_k4(1) + X_k4(12)*F1t_l_k4(1) + X_k4(20)*v_orbit1_k4(1);
    v1y_k4 = X_k4(11)*F1n_l_k4(2) + X_k4(12)*F1t_l_k4(2) + X_k4(20)*v_orbit1_k4(2);
    v1z_k4 = X_k4(11)*F1n_l_k4(3) + X_k4(12)*F1t_l_k4(3) + X_k4(20)*v_orbit1_k4(3);
    v2x_k4 = X_k4(13)*F2n_l_k4(1) + X_k4(14)*F2t_l_k4(1) + X_k4(21)*v_orbit2_k4(1);
    v2y_k4 = X_k4(13)*F2n_l_k4(2) + X_k4(14)*F2t_l_k4(2) + X_k4(21)*v_orbit2_k4(2);
    v2z_k4 = X_k4(13)*F2n_l_k4(3) + X_k4(14)*F2t_l_k4(3) + X_k4(21)*v_orbit2_k4(3);
    
    F1_k4=X_k4(1)*F1n_l_k4 + F1t*F1t_l_k4 + X_k4(19)*v_orbit1_k4;
    F2_k4=X_k4(2)*F2n_l_k4 + F2t*F2t_l_k4 + X_k4(20)*v_orbit2_k4;
    F0_k4 = [X_k4(8),X_k4(9),X_k4(10)]';
    dphi1_k4 = X_k4(6);
    dphi2_k4 = X_k4(7);
    w0x_k4 = X_k4(3);
    w0y_k4 = X_k4(4);
    w0z_k4 = X_k4(5);

%% update velocity
    v0x_k = (v0x_k1 + 2*v0x_k2 + 2*v0x_k3 + v0x_k4)/6;
    v0y_k = (v0y_k1 + 2*v0y_k2 + 2*v0y_k3 + v0y_k4)/6;
    v0z_k = (v0z_k1 + 2*v0z_k2 + 2*v0z_k3 + v0z_k4)/6;
    v1x_k = (v1x_k1 + 2*v1x_k2 + 2*v1x_k3 + v1x_k4)/6;
    v1y_k = (v1y_k1 + 2*v1y_k2 + 2*v1y_k3 + v1y_k4)/6;
    v1z_k = (v1z_k1 + 2*v1z_k2 + 2*v1z_k3 + v1z_k4)/6;
    v2x_k = (v2x_k1 + 2*v2x_k2 + 2*v2x_k3 + v2x_k4)/6;
    v2y_k = (v2y_k1 + 2*v2y_k2 + 2*v2y_k3 + v2y_k4)/6;
    v2z_k = (v2z_k1 + 2*v2z_k2 + 2*v2z_k3 + v2z_k4)/6;
    dphi1_k = (dphi1_k1 + 2*dphi1_k2 + 2*dphi1_k3 + dphi1_k4)/6;
    dphi2_k = (dphi2_k1 + 2*dphi2_k2 + 2*dphi2_k3 + dphi2_k4)/6;
    w0x_k = (w0x_k1 + 2*w0x_k2 + 2*w0x_k3 + w0x_k4)/6;
    w0y_k = (w0y_k1 + 2*w0y_k2 + 2*w0y_k3 + w0y_k4)/6;
    w0z_k = (w0z_k1 + 2*w0z_k2 + 2*w0z_k3 + w0z_k4)/6;
    F1_k = (F1_k1 + 2*F1_k2 + 2*F1_k3 + F1_k4)/6;
    F2_k = (F2_k1 + 2*F2_k2 + 2*F2_k3 + F2_k4)/6;
    F0_k = (F0_k1 + 2*F0_k2 + 2*F0_k3 + F0_k4)/6;
    Lw =  [w0x_k; w0y_k; w0z_k];
    Ldr = dr0_l;
    Ldy = dyaw_l;
    Ldp = dp_l;

    L_M = [Ldr(1), Ldy(1), Ldp(1);  Ldr(2), Ldy(2), Ldp(2);  Ldr(3), Ldy(3), Ldp(3)];
    Bw = L_M \ Lw;

    w0x_k=Bw(1);
    w0y_k=Bw(2);
    w0z_k=Bw(3);

    
%% update position
% 1: By velocity
    r0_l(1)=r0_l(1)+dt*v0x_k;
    r0_l(2)=r0_l(2)+dt*v0y_k;
    r0_l(3)=r0_l(3)+dt*v0z_k;
%     r1_l(1)=r1_l(1)+dt*v1x_k;
%     r1_l(2)=r1_l(2)+dt*v1y_k;
%     r1_l(3)=r1_l(3)+dt*v1z_k;
%     r2_l(1)=r2_l(1)+dt*v2x_k;
%     r2_l(2)=r2_l(2)+dt*v2y_k;
%     r2_l(3)=r2_l(3)+dt*v2z_k;

% 2: By model dynamics
    phi1=phi1+dt*dphi1_k;
    phi2=phi2+dt*dphi2_k;

    w0x_k = +w0x_k;
    w0y_k = +w0y_k;
    w0z_k = +w0z_k;
%     abwn(:,i)=[w0x_k; w0y_k; w0z_k];

    dr0_l            =  rotate(dr0_l,dr0_l,w0x_k*dt);
    deye_l           =  rotate(deye_l,dr0_l,w0x_k*dt);
    dyaw_l           =  rotate(dyaw_l,dr0_l,w0x_k*dt);
    dp_l             =  rotate(dp_l,dr0_l,w0x_k*dt);
    v_orbit1         =  rotate(v_orbit1,dr0_l,w0x_k*dt);
    v_orbit2         =  rotate(v_orbit2,dr0_l,w0x_k*dt);
    vn1_0            =  rotate(vn1_0,dr0_l,w0x_k*dt);
    vn2_0            =  rotate(vn2_0,dr0_l,w0x_k*dt);
    ex_b             =  rotate(ex_b,dr0_l,w0x_k*dt);
    ey_b             =  rotate(ey_b,dr0_l,w0x_k*dt);
    ez_b             =  rotate(ez_b,dr0_l,w0x_k*dt);

    dr0_l            =  rotate(dr0_l,dyaw_l,w0y_k*dt);
    deye_l           =  rotate(deye_l,dyaw_l,w0y_k*dt);
    dyaw_l           =  rotate(dyaw_l,dyaw_l,w0y_k*dt);
    dp_l             =  rotate(dp_l,dyaw_l,w0y_k*dt);
    v_orbit1         =  rotate(v_orbit1,dyaw_l,w0y_k*dt);
    v_orbit2         =  rotate(v_orbit2,dyaw_l,w0y_k*dt);
    vn1_0            =  rotate(vn1_0,dyaw_l,w0y_k*dt);
    vn2_0            =  rotate(vn2_0,dyaw_l,w0y_k*dt);
    ex_b             =  rotate(ex_b,dyaw_l,w0y_k*dt);
    ey_b             =  rotate(ey_b,dyaw_l,w0y_k*dt);
    ez_b             =  rotate(ez_b,dyaw_l,w0y_k*dt);

    dr0_l            =  rotate(dr0_l,dp_l,w0z_k*dt);
    deye_l           =  rotate(deye_l,dp_l,w0z_k*dt);
    dyaw_l           =  rotate(dyaw_l,dp_l,w0z_k*dt);
    dp_l             =  rotate(dp_l,dp_l,w0z_k*dt);
    v_orbit1         =  rotate(v_orbit1,dp_l,w0z_k*dt);
    v_orbit2         =  rotate(v_orbit2,dp_l,w0z_k*dt);
    vn1_0            =  rotate(vn1_0,dp_l,w0z_k*dt);
    vn2_0            =  rotate(vn2_0,dp_l,w0z_k*dt);
    ex_b             =  rotate(ex_b,dp_l,w0z_k*dt);
    ey_b             =  rotate(ey_b,dp_l,w0z_k*dt);
    ez_b             =  rotate(ez_b,dp_l,w0z_k*dt);

    F1n_l            =  rotate(vn1_0,v_orbit1,phi1);
    F2n_l            =  rotate(vn2_0,v_orbit2,phi2);
    F1t_l            =  rotate(F1n_l,v_orbit1,pi/2);
    F2t_l            =  rotate(F2n_l,v_orbit2,pi/2);

    phi11=atan2(sin(phi1)*Ra1,cos(phi1)*Rb1);
    phi22=atan2(sin(phi2)*Ra2,cos(phi2)*Rb2);
    R1=sqrt(Ra1*Ra1*cos(phi11)*cos(phi11) + Rb1*Rb1*sin(phi11)*sin(phi11));
    R2=sqrt(Ra2*Ra2*cos(phi22)*cos(phi22) + Rb2*Rb2*sin(phi22)*sin(phi22));

    r01=-l1*ex_b + h1*ez_b;
    r02= l2*ex_b + h2*ez_b;

    r1_l=r0_l+r01+R1*F1n_l;
    r2_l=r0_l+r02+R2*F2n_l;

    r1_lp=r0_l+r01/b1+R1/b2*F1n_l;
    r2_lp=r0_l+r02/b1+R2/b2*F2n_l;
    
    %Calculate work done
    %1:
    torq_old = 8*pi*mu*a0*a0*a0*Lw;
    V_0dt = [v0x_k,v0y_k,v0z_k]'*dt;
    w_0dt = Lw*dt;
    W1 = dot(-F0_k,V_0dt);
    W2 = dot(torq_old,w_0dt);
    W_all = W1+W2;
    Work = Work+W_all;
    w1all = w1all+W1;
    w2all = w2all+W2;
    W_fla = dot(F1_k,[v1x_k,v1y_k,v1z_k]'*dt) + dot(F2_k,[v2x_k,v2y_k,v2z_k]'*dt);
    W_motor = W_motor+W_fla;
    wfla = [wfla,W_fla];
    wleft = [wleft,dot(F1_k,[v1x_k,v1y_k,v1z_k]'*dt)];
    wright = [wright,dot(F2_k,[v2x_k,v2y_k,v2z_k]'*dt)];
    dis_e = norm(r0_old-r0_l);
    d_ff=[d_ff,dis_e];
    V0=[V0;[v0x_k,v0y_k,v0z_k]];
    d_vector = (r0_old-r0_l);
    ang_d1 = atan2(norm(cross(F1_k,d_vector)), dot(F1_k,d_vector));
    ang_d2 = atan2(norm(cross(F2_k,d_vector)), dot(F2_k,d_vector));
    ang1 = [ang1,ang_d1];
    ang2 = [ang2,ang_d2];
end

info = {ar0n,aden,adr0n,aphi1n,aphi2n,ang_d1,ang_d2};
end


function [vector1]=rotate(vector,axis,theta)
k_hat = [0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0 ];
R = eye(3)+ k_hat.*sin(theta) + k_hat*k_hat.*(1-cos(theta));
vector1 = R*vector;
end