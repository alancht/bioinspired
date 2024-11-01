clear 
clc
close all
mu=(10.6/10)^2;
b1 = 0.8;
b2 = 1.5;
afact=0.5;
minV = round(0.42,2);
maxV = round(0.7,2);

incre = 0.01;
rand_init = (maxV-minV).*rand(2,1) + minV;
l1p = round(rand_init(1),2);
l2p = round(rand_init(2),2);
h1p = sqrt(0.5-l1p^2);
h2p = sqrt(0.5-l2p^2);
%

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


Ra1=0.1572*b2;
Rb1=0.2803*b2;
Ra2=0.1649*b2;
Rb2=0.2906*b2;


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
nt = 2000;
dt = 0.05;
r1ns=R1ns;
r2ns=R2ns;




state_num = int32(((maxV-minV)*(1/incre)+1)^2);

action1 = [-0.05,-0.03,-0.01,0,0.01,0.03,0.05];

action2 = action1;
action_all = [];
for i=1:length(action1)
    temp1 = [action1(i)];
    for j=1:length(action2)
        temp = [temp1,action2(j)];
        action_all = [action_all;temp];
        
    end
end

[action_num,pass1] = size(action_all);


l1s = round(minV:incre:maxV,2);
l2s = round(minV:incre:maxV,2);
numb = 1:length(l1s)^2;
Nstate = reshape(numb,[length(l1s),length(l1s)]);
[idxa,idxb] = find(Nstate==5);


eps_count = 0;
start_eps = 500;
eps_limit = 3500;
epsilon_1 = 0.2;
epsilon_2 = 0.01;
alpha=1;
gamma=0.98;

check_eps = 500;
checkpoint = 500;
wmin = 1000;
wmax = -1000;
dmin = 1000;
dmax = -1000;


parfor alph=0:10
        eps_r = [];
        Q=-1000*ones(state_num,action_num);
        R1ns=r1ns;
        R2ns=r2ns;
    for eps_count=1:(eps_limit+check_eps)
        afact=alph/10;
        a_list = [];
        s_list = [];
        xn_l = [];
        wlist = [];
        dlist = [];
      
    
        rand_init = (maxV-minV).*rand(2,1) + minV;
        l1p = round(rand_init(1),2);
        l2p = round(rand_init(2),2);

        h1p = sqrt(0.5-l1p^2);
        h2p = sqrt(0.5-l2p^2);
        
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
        action=ceil(rand*action_num);
    
        %Find S
        idX = int32(find(l1s==l1p));
        idY = int32(find(l2s==l2p));
        s = Nstate(idX,idY);
    
        reward_all = [];

        if eps_count>eps_limit 
            epsilon = 0;
        else
            epsilon = (epsilon_2-epsilon_1)/eps_limit*eps_count+epsilon_1;
        end
        tic
        for i=1:100
            bound = 0;
            randcheck=rand;
            if randcheck>epsilon
                action=find(Q(s,:)==max(Q(s,:)));  
    
            end
            if (length(action)>1 || randcheck<=epsilon) 
                randidx = rand;
                action=ceil(action_num*randidx);
            end

            a_list = [a_list,action];
            step1 = action_all(action,1);
            step2 = action_all(action,2);
    
            l1t = round(l1p+step1,2);
            l2t = round(l2p+step2,2);
            if l1t<minV || l1t>maxV || l2t<minV || l2t>maxV
                bound = 1;
            end
    
            l1p = round(clip(l1p+step1,minV,maxV),2);
            l2p = round(clip(l2p+step2,minV,maxV),2);
            h1p = sqrt(0.5-l1p^2);
            h2p = sqrt(0.5-l2p^2);
            

            idX = int32(find(l1s==l1p));
            idY = int32(find(l2s==l2p));
            s_next = Nstate(idX,idY);
    
            r0_old = r0_l;
            temp = [l1p,l2p];
            s_list = [s_list;temp];
            [r0_l,phi1,phi2,dr0_l,deye_l,dyaw_l,dp_l,ex_b,ey_b,ez_b,R1ns,R2ns,W0,w1all,w2all,W_motor,info,dff,V0]...
                =swimstep(r0_l,phi1,phi2,dr0_l,deye_l,dyaw_l,dp_l,ex_b,ey_b,ez_b,R1ns,R2ns,l1p,l2p,h1p,h2p,nt,dt);
            d_all = sum(dff);
            dP = norm(r0_l-r0_old);
            
            Punit=(r0_l-r0_old)/norm(r0_l-r0_old);
            V_parallel=V0*Punit.*Punit';
            Vpdt=dt.*V_parallel;
            F_Vis=6*pi*mu*0.5.*V0;
            FVdt=F_Vis.*Vpdt;
            work_eff=sum(sum(FVdt));

        
    

            
            W_stad = -work_eff;
            D_stad = dP;
            reward = afact*W_stad+(1-afact)*D_stad;
            if bound == 1
                reward = -2*abs(reward);
            end
            Q(s,action)=Q(s,action)+alpha*(reward+gamma*max(Q(s_next,:))-Q(s,action));
            reward_all = [reward_all,reward];

            s = s_next;
           
            wlist = [wlist,work_eff];
            dlist = [dlist,dP];

    
        end
        cul_reward = sum(reward_all);
    
    
        disp(s_list(end,:))
        eps_r = [eps_r,cul_reward];
        toc
        
    
    end
end


function clip_res = clip(input, minVal, maxVal)
    input(input < minVal) = minVal;
    input(input > maxVal) = maxVal;
    clip_res = input;
end

