clear 
clc
close all

mu=(10.6/10)^2;
b1 = 0.8;
b2 = 1.5;


%

r0_l=[0;0;0]; 
ex_b=[1;0;0];          
ey_b=[0;1;0];
ez_b=[0;0;1];
turn_d = [0;1;0];
phi1=0*pi;            
phi2=0*pi;   

dr0_l=[0 ; 0;  1];                      
deye_l=[sin(pi/4);  -sin(pi/4);  0];    
dyaw_l=deye_l;                          
dp_l =[-sin(pi/4);  -sin(pi/4);  0];    


Ra1=0.1572*b2;
Rb1=0.2803*b2;
Ra2=0.1649*b2;
Rb2=0.2906*b2;

%% dr/dtheta
R1n=zeros(1,360);
R2n=zeros(1,360);
R1ns=zeros(1,360);
R2ns=zeros(1,360);
for i = 0:pi/180:2*pi
    phi11=atan2(sin(i)*Ra1,cos(i)*Rb1);
    phi22=atan2(sin(i)*Ra2,cos(i)*Rb2);
    i1=round(i/(pi/180)+1);
    R1n(1,i1)=sqrt(Ra1*Ra1*cos(phi11)*cos(phi11) + Rb1*Rb1*sin(phi11)*sin(phi11));
    R2n(1,i1)=sqrt(Ra2*Ra2*cos(phi22)*cos(phi22) + Rb2*Rb2*sin(phi22)*sin(phi22));
end
for i = 1:360
    R1ns(1,i)=(R1n(1,i+1)-R1n(1,i))*180/pi;
    R2ns(1,i)=(R2n(1,i+1)-R2n(1,i))*180/pi;
end
nt = 1500;
dt = 0.05;
r1ns=R1ns;
r2ns=R2ns;

eps_count = 0;

eps_limit = 60;

epsilon = 0.02;
alpha=1;
gamma=0.95;

check_eps = 0;

turn_d = [1;0;0];

eps_r = [];
state_num = 21*4;
action_num = 3;
action_all = [0.15,0;0,0.15;0,0];
Q=rand(state_num,action_num);
R1ns=r1ns;
R2ns=r2ns;
pre_rew = -100;
for eps_count=1:(eps_limit+check_eps)

    a_list = [];
    s_list = [];
    xn_l = [];




    l1i = 0.45;
    l2i = 0.42;
    h1p = sqrt(0.5-l1i^2);
    h2p = sqrt(0.5-l2i^2);
    new_exb = rotate(ex_b,ez_b,pi/4);
    new_eyb = rotate(ey_b,ez_b,pi/4);

    vec = dot(ez_b,turn_d);
    if dot(new_exb,turn_d)>=0 && dot(new_eyb,turn_d)>=0
        gs=1;
    elseif dot(new_exb,turn_d)>=0 && dot(new_eyb,turn_d)<0
        gs=2;
    elseif dot(new_exb,turn_d)<0 && dot(new_eyb,turn_d)<0
        gs=3;
    elseif dot(new_exb,turn_d)<0 && dot(new_eyb,turn_d)>=0
        gs=4;
    end
    pointer = floor((vec+1)*10)+1;
    s = (gs-1)*21+pointer;
    r0_l=[0;0;0]; 
    ex_b=[1;0;0];         
    ey_b=[0;1;0];
    ez_b=[0;0;1];
    
    phi1=0*pi;            
    phi2=0*pi;   
    
    dr0_l=[0 ; 0;  1];                      
    deye_l=[sin(pi/4);  -sin(pi/4);  0];    
    dyaw_l=deye_l;                          
    dp_l =[-sin(pi/4);  -sin(pi/4);  0];    


   

    reward_all = [];
    S_rec=[];
    A_rec=[];


    tic
    for i=1:30
        bound = 0;
        randcheck=rand;
        if randcheck>epsilon
            action=find(Q(s,:)==max(Q(s,:)));  

        end
        if (length(action)>1 || randcheck<=epsilon) 
            randidx = rand;
            action=ceil(action_num*randidx);
        end
        S_rec=[S_rec;s];
        A_rec=[A_rec;action];

        step1 = action_all(action,1);
        step2 = action_all(action,2);

        l1p = round(l1i+step1,2);
        l2p = round(l2i+step2,2);
        h1p = sqrt(0.5-l1p^2);
        h2p = sqrt(0.5-l2p^2);
        

        new_exb = rotate(ex_b,ez_b,pi/4);
        new_eyb = rotate(ey_b,ez_b,pi/4);

        vec = dot(ez_b,turn_d);
        if dot(new_exb,turn_d)>=0 && dot(new_eyb,turn_d)>=0
            gs=1;
        elseif dot(new_exb,turn_d)>=0 && dot(new_eyb,turn_d)<0
            gs=2;
        elseif dot(new_exb,turn_d)<0 && dot(new_eyb,turn_d)<0
            gs=3;
        elseif dot(new_exb,turn_d)<0 && dot(new_eyb,turn_d)>=0
            gs=4;
        end
        pointer = floor((vec+1)*10)+1;
        s_next = (gs-1)*21+pointer;

        s_list = [s_list,s_next];
        a_list = [a_list,action];
        r0_old = r0_l;
 
        [r0_l,phi1,phi2,dr0_l,deye_l,dyaw_l,dp_l,ex_b,ey_b,ez_b,R1ns,R2ns,W0,w1all,w2all,W_motor,info,dff]...
            =swimstep(r0_l,phi1,phi2,dr0_l,deye_l,dyaw_l,dp_l,ex_b,ey_b,ez_b,R1ns,R2ns,l1p,l2p,h1p,h2p,nt,dt);
        d_all = sum(dff);
        dP = norm(r0_l-r0_old);
        reward = dot(ez_b,turn_d);
    

        %Update Q matrix
        

        reward_all = [reward_all,reward];

        s = s_next;



    end
    i=i+1;
    Rew=0;
    cul_reward = sum(reward_all);
    for Tt=length(reward_all):-1:1
        Rew=Rew*gamma+reward_all(Tt);
        up_s=S_rec(Tt);
        up_a=A_rec(Tt);

        
        if abs(Q(up_s,up_a))<abs(Rew)
            Q(up_s,up_a)=Rew;
        end

    end
    if pre_rew<cul_reward
        pre_rew = cul_reward;
    end
    





    toc

   

end


function [vector1]=rotate(vector,axis,theta)
k_hat = [0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0 ];
R = eye(3)+ k_hat.*sin(theta) + k_hat*k_hat.*(1-cos(theta));
vector1 = R*vector;
end