%%%%%%%%%%%% Optimal control with delay %%%%%%%%%%%%%
function [z1,z2,z3,z4,z5,S,E,I,H,R,u11,v1,w1,z11,tt] =  compute_optimal_control_delay(taux1,taux2,taux3,taux4)
% NN # of days

% clc;
%clear all;
global eta beta mu N alpha theta delta1 lamda gamma delta2

eta=2000; mu=0.02; delta1=0.2; delta2=0.03; beta=0.6; alpha=0.4; %beta=0.045
theta=0.2; lamda=0.01; gamma=0.1;
A=10^(+7); B=10^6; C=5*10^6; D=10^6;
S0=2*10^9; I0=10^5; E0=10^6; H0=0; R0=0;

delta=0.001;

S=S0; E=E0; I=I0; H=H0; R=R0;
MaxTime=100;
test=-1;
NN=100;
N=S+E+I+H+R;


MaxTime = MaxTime+taux4; %REVOIRE APRES URGENT !!!

NN=100+taux4; % pour prendre en compte retard

t=linspace(0,MaxTime,NN+1); % linspace(0,600,101) % get 101 points (pas =6)


%  I is the interval [0,N+1]
h=MaxTime/NN; %REVOIRE APRES URGENT !!!

h2=h/2;

u=zeros(1,NN+1+2*taux4);
v=zeros(1,NN+1+2*taux4);
w=zeros(1,NN+1+2*taux4);
z=zeros(1,NN+1+2*taux4);

%% initialize variables without control
z1=zeros(1,NN+1); z2=zeros(1,NN+1); z3=zeros(1,NN+1);z4=zeros(1,NN+1);z5=zeros(1,NN+1);
z1(:,1:taux4+1)=S0; z2(:,1:taux4+1)=E0; z3(:,1:taux4+1)=I0; z4(:,1:taux4+1)=H0; z5(:,1:taux4+1)=R0;

%% initialize variables with control
S=zeros(1,NN+1+2*taux4); E=zeros(1,NN+1+2*taux4); I=zeros(1,NN+1+2*taux4); H=zeros(1,NN+1+2*taux4);
% declare SEIR variables considering delay taux1
S_taux1=zeros(1,NN+1+2*taux4); E_taux1=zeros(1,NN+1+2*taux4);
% declare SEIR variables considering delay taux2
I_taux2=zeros(1,NN+1+2*taux4);
% declare SEIR variables considering delay taux3
E_taux3=zeros(1,NN+1+2*taux4);
% declare SEIR variables considering delay taux4
H_taux4=zeros(1,NN+1+2*taux4);

R=zeros(1,NN+1+2*taux4);

S(:,1:taux4+1)=S0; E(:,1:taux4+1)=E0; I(:,1:taux4+1)=I0; H(:,1:taux4+1)=H0; R(:,1:taux4+1)=R0;


% initialize sigma variables
sigma1=zeros(1,NN+1+2*taux4); sigma2=zeros(1,NN+1+2*taux4); sigma3=zeros(1,NN+1+2*taux4);
sigma4=zeros(1,NN+1+2*taux4); sigma5=zeros(1,NN+1+2*taux4);
% initialize sigma variables considering delay taux
sigma1_taux1=zeros(1,NN+1+2*taux4);
sigma4_taux2=zeros(1,NN+1+2*taux4);
sigma4_taux3=zeros(1,NN+1+2*taux4);
sigma5_taux4=zeros(1,NN+1+2*taux4);


%condition de transversalité A REVOIR
sigma2(NN+taux4+1)=1; sigma3(NN+taux4+1)=1; sigma4(NN+taux4+1)=1;

x=[S E I H R]'; Z=[z1 z2 z3 z4 z5]';

sigma=[sigma1 sigma2 sigma3 sigma4 sigma5]';
u1=u';

tt = 1:1:NN+2*taux4+1;
%tt=linspace(0,MaxTime,NN+1);

for i=1+taux4:NN+2*taux4 % loop for NN+taux4 times
    %  for i=1:NN % loop for NN+taux4 times
    k11=eta - beta*z1(i)*z2(i)/N - mu*z1(i) ;
    k12=beta*z1(i)*z2(i)/N - (mu+alpha+theta)*z2(i) ;
    k13=alpha*z2(i) - (mu+lamda+delta1)*z3(i) ;
    k14=lamda*z3(i) - (mu+gamma+delta2)*z4(i) ;
    k15=gamma*z4(i) + theta*z2(i) - mu*z5(i) ;
    
    k21=eta - beta*(z1(i)+h2*k11)*(z2(i)+h2*k12)/N - mu*(z1(i)+h2*k11) ;
    k22=beta*(z1(i)+h2*k11)*(z2(i)+h2*k12)/N - (mu+alpha+theta)*(z2(i)+h2*k12) ;
    k23=alpha*(z2(i)+h2*k12) - (mu+lamda+delta1)*(z3(i)+h2*k13) ;
    k24=lamda*(z3(i)+h2*k13) - (mu+gamma+delta2)*(z4(i)+h2*k14) ;
    k25=gamma*(z4(i)+h2*k14) + theta*(z2(i)+h2*k12) - mu*(z5(i)+h2*k15) ;
    
    k31=eta - beta*(z1(i)+h2*k21)*(z2(i)+h2*k22)/N - mu*(z1(i)+h2*k21) ;
    k32=beta*(z1(i)+h2*k21)*(z2(i)+h2*k22)/N - (mu+alpha+theta)*(z2(i)+h2*k22) ;
    k33=alpha*(z2(i)+h2*k22) - (mu+lamda+delta1)*(z3(i)+h2*k23) ;
    k34=lamda*(z3(i)+h2*k23) - (mu+gamma+delta2)*(z4(i)+h2*k24) ;
    k35=gamma*(z4(i)+h2*k24) + theta*(z2(i)+h2*k22) - mu*(z5(i)+h2*k25) ;
    
    k41=eta - beta*(z1(i)+h*k31)*(z2(i)+h*k32)/N - mu*(z1(i)+h*k31) ;
    k42=beta*(z1(i)+h*k31)*(z2(i)+h*k32)/N - (mu+alpha+theta)*(z2(i)+h*k32) ;
    k43=alpha*(z2(i)+h*k32) - (mu+lamda+delta1)*(z3(i)+h*k33) ;
    k44=lamda*(z3(i)+h*k33) - (mu+gamma+delta2)*(z4(i)+h*k34) ;
    k45=gamma*(z4(i)+h*k34) + theta*(z2(i)+h*k32) - mu*(z5(i)+h*k35) ;
    
    z1(i+1)=z1(i)+(h/6)*(k11+2*k21+2*k31+k41);
    z2(i+1)=z2(i)+(h/6)*(k12+2*k22+2*k32+k42);
    z3(i+1)=z3(i)+(h/6)*(k13+2*k23+2*k33+k43);
    z4(i+1)=z4(i)+(h/6)*(k14+2*k24+2*k34+k44);
    z5(i+1)=z5(i)+(h/6)*(k15+2*k25+2*k35+k45);
end

%% WITH CONTROL
count = 1;
%% looking for solutions with control

while(test<0 && count<=20)
    
    oldu=u;
    oldu1=u1; oldv=v; oldw=w;  oldz=z;
    oldx=x;
    oldsigma=sigma;
    
    for i=1+taux4:NN+taux4   % to obtain S,E,I,and R with dimension [(NN+taux4)x1]
        k11=eta - mu*S(i) - beta*S(i)*E(i)/N + beta*(u(i-taux1))*S(i-taux1)*E(i-taux1)/N   ;
        k12=beta*(1-u(i))*S(i)*E(i)/N - (mu+alpha+theta)*E(i) - w(i)*E(i) ;
        k13=alpha*E(i) - (mu+lamda+delta1)*I(i) - v(i)*I(i) ;
        k14=lamda*I(i) - (mu+gamma+delta2)*H(i) +w(i-taux3)*E(i-taux3) + v(i-taux2)*I(i-taux2) -z(i)*H(i) ;
        k15=gamma*H(i) + theta*E(i) - mu*R(i) + z(i-taux4)* H(i-taux4);
        
        
        k21=eta - mu*(S(i)+h2*k11) - beta*(S(i)+h2*k11)*(E(i)+h2*k12)/N + beta*0.5*(u(i-taux1)+u(i-taux1+1))*(S(i-taux1)+h2*k11)*(E(i-taux1)+h2*k12)/N;
        
        k22=beta*(1-0.5*(u(i)+u(i+1))).*(S(i)+h2*k11)*(E(i)+h2*k12)/N - (mu+alpha+theta)*(E(i)+h2*k12) - 0.5*(w(i)+w(i+1))*(E(i)+h2*k12) ;
        
        k23=alpha*(E(i)+h2*k12) - (mu+lamda+delta1)*(I(i)+h2*k13) - 0.5*(v(i)+v(i+1))*(I(i)+h2*k13) ;
        
        k24=lamda*(I(i)+h2*k13) - (mu+gamma+delta2)*(H(i)+h2*k14) + 0.5*(v(i-taux2)+v(i-taux2+1))*(I(i-taux2)+h2*k13) + 0.5*(w(i-taux3)+w(i-taux3+1))*(E(i-taux3)+h2*k12) - 0.5*(z(i)+z(i+1))*(H(i)+h2*k14) ;
        k25=gamma*(H(i)+h2*k14) + theta*(E(i)+h2*k12) - mu*(R(i)+h2*k15) + 0.5*(z(i-taux4)+z(i-taux4+1))*(H(i-taux4)+h2*k14) ;
        
        k31=eta - mu*(S(i)+h2*k21) - beta*(S(i)+h2*k21)*(E(i)+h2*k22)/N + beta*0.5*(u(i-taux1)+u(i-taux1+1))*(S(i-taux1)+h2*k21)*(E(i-taux1)+h2*k22)/N;
        k32=beta*(1-0.5*(u(i)+u(i+1)))*(S(i)+h2*k21)*(E(i)+h2*k22)/N - (mu+alpha+theta)*(E(i)+h2*k22) - 0.5*(w(i)+w(i+1))*(E(i)+h2*k22) ;
        k33=alpha*(E(i)+h2*k22) - (mu+lamda+delta1)*(I(i)+h2*k23) - 0.5*(v(i)+v(i+1))*(I(i)+h2*k23) ;
        k34=lamda*(I(i)+h2*k23) - (mu+gamma+delta2)*(H(i)+h2*k24) + 0.5*(v(i-taux2)+v(i-taux2+1))*(I(i-taux2)+h2*k23) + 0.5*(w(i-taux3)+w(i-taux3+1))*(E(i-taux3)+h2*k22) - 0.5*(z(i)+z(i+1))*(H(i)+h2*k24) ;
        k35=gamma*(H(i)+h2*k24) + theta*(E(i)+h2*k22) - mu*(R(i)+h2*k25) + 0.5*(z(i-taux4)+z(i-taux4+1))*(H(i-taux4)+h2*k24) ;
        
        k41=eta - mu*(S(i)+h2*k31) - beta*(S(i)+h2*k31)*(E(i)+h2*k32)/N + beta*u(i-taux1+1)*(S(i-taux1)+h2*k31)*(E(i-taux1)+h2*k32)/N;
        k42=beta*(1-u(i+1))*(S(i)+h2*k31)*(E(i)+h2*k32)/N - (mu+alpha+theta)*(E(i)+h2*k32) - w(i+1)*(E(i)+h2*k32) ;
        k43=alpha*(E(i)+h2*k32) - (mu+lamda+delta1)*(I(i)+h2*k33) - v(i+1)*(I(i)+h2*k33) ;
        k44=lamda*(I(i)+h2*k33) - (mu+gamma+delta2)*(H(i)+h2*k34) + v(i-taux2+1)*(I(i-taux2)+h2*k33) + w(i-taux3+1)*(E(i-taux3)+h2*k32) - z(i+1)*(H(i)+h2*k34) ;
        k45=gamma*(H(i)+h2*k34) + theta*(E(i)+h2*k32) - mu*(R(i)+h2*k35) + z(i-taux4+1)*(H(i-taux4)+h2*k34) ;
        
        %j= i-taux4;
        S(i+1)=S(i)+(h/6)*(k11+2*k21+2*k31+k41);% [(NN+taux4)x1]
        E(i+1)=E(i)+(h/6)*(k12+2*k22+2*k32+k42);% [(NN+taux4)x1]
        I(i+1)=I(i)+(h/6)*(k13+2*k23+2*k33+k43);% [(NN+taux4)x1]
        H(i+1)=H(i)+(h/6)*(k14+2*k24+2*k34+k44);% [(NN+taux4)x1]
        R(i+1)=R(i)+(h/6)*(k15+2*k25+2*k35+k45);% [(NN+taux4)x1]
        
        
    end
    
    fprintf('counting iterations %f \n',count);
    
    % SEIHR with delays taux
    S_taux1(1:end-taux1) = S(1+taux1 : end); % [(NN+taux4-taux1)x1]
    I_taux2(1:end-taux2) =  I(1+taux2 : end); % [(NN+taux4-taux2)x1]
    E_taux1(1:end-taux1) = E(1+taux1: end); % [(NN+taux4-taux4)x1]
    E_taux3(1:end-taux3) = E(1+taux3: end); % [(NN+taux4-taux3)x1]
    H_taux4(1:end-taux4)= H(1+taux4 : end); %[(NN)x1]
    
    for i=1+taux4:NN+taux4
        j=NN+2*taux4+2-i;
        
        k11=(mu+beta*E(j)/N)*sigma1(j)-(beta*(1-u(j))*E(j)/N)*sigma2(j) - beta*double(j<=(NN+1-taux1) && j>=0) * ...
            sigma1(j+taux1)*u(j+taux1)*E(j+taux1)/N;% [(NN+taux4)x1]
        k12=-1+(beta*S(j)/N)*sigma1(j)-(beta*(1-u(j))*S(j)/N -mu-alpha-theta -w(j)) ...
            *sigma2(j) - alpha*sigma3(j) - theta*sigma5(j) - double(j<=(NN+1-taux3) && j>=0) * ...
            sigma4(j+taux3)*w(j+taux3);% [(NN+taux4)x1] % + beta *u*s/N*sigma1 is deleted
        k13=-1+(mu+lamda+delta1+v(j))*sigma3(j) - lamda*sigma4(j) - double(j<=(NN+1-taux2) && j>=0) * ...
            sigma4(j+taux2)*v(j+taux2);% [(NN+taux4)x1]
        k14=-1+(mu+gamma+delta2+z(j))*sigma4(j) - gamma*sigma5(j) - double(j<(NN-taux4) && j>0) * ...
            sigma5(j+taux4)*z(j+taux4);% [(NN+taux4)x1]
        k15=mu*sigma5(j);% [(NN+taux4)x1]
        
        
        k21=(mu+beta*0.5*(E(j)+E(j-1))/N)*(sigma1(j)-h2*k11)-...
            (beta*(1-0.5*(u(j)+u(j-1)))*0.5*(E(j)+E(j-1))/N)*(sigma2(j)-h2*k12) - beta* ...
            double(j<=(NN+1-taux1) && j>=0) * ...
            (sigma1(j+taux1)-h2*k11)*0.5*(u(j+taux1)+u(j+taux1-1))*0.5*(E(j+taux1)+E(j+taux1-1))/N;
        
        
        k22=-1+(beta*0.5*(S(j)+S(j-1))/N)*(sigma1(j)-h2*k11)-...
            (beta*(1-0.5*(u(j)+u(j-1)))*0.5*(S(j)+S(j-1))/N-mu-alpha-theta-...
            0.5*(w(j)+w(j-1)))*(sigma2(j)-h2*k12) - alpha*(sigma3(j)-h2*k13) -...
            theta*(sigma5(j)-h2*k15) - double(j<=(NN+1-taux3) && j>=0) * ...
            (sigma4(j+taux3)-h2*k14)*0.5*(w(j+taux3)+w(j+taux3-1));
        
        k23=-1+(mu+lamda+delta1+0.5*(v(j)+v(j-1)))*(sigma3(j)-h2*k13) - ...
            lamda*(sigma4(j)-h2*k14) - double(j<=(NN+1-taux2) && j>=0) * ...
            (sigma4(j+taux2)-h2*k14)*0.5*(v(j+taux2)+v(j+taux2-1));
        k24=-1+(mu+gamma+delta2+0.5*(z(j)+z(j-1)))*(sigma4(j)-h2*k14) - ...
            gamma*(sigma5(j)-h2*k15) - double(j<(NN-taux4) && j>0) * ...
            (sigma5(j+taux4)-h2*k15)*0.5*(z(j+taux4)+z(j+taux4-1));
        k25=mu*(sigma5(j)-h2*k15);
        
        
        k31=(mu+beta*0.5*(E(j)+E(j-1))/N)*(sigma1(j)-h2*k21)-...
            (beta*(1-0.5*(u(j)+u(j-1)))*0.5*(E(j)+E(j-1))/N)*(sigma2(j)-h2*k22) - beta* ...
            double(j<=(NN+1-taux1) && j>=0) * ...
            (sigma1(j+taux1)-h2*k21)*0.5*(u(j+taux1)+u(j+taux1-1))*0.5*(E(j+taux1)+E(j+taux1-1))/N;
        
        k32=-1+(beta*0.5*(S(j)+S(j-1))/N)*(sigma1(j)-h2*k21)-...
            (beta*(1-0.5*(u(j)+u(j-1)))*0.5*(S(j)+S(j-1))/N-mu-alpha-theta-...
            0.5*(w(j)+w(j-1)))*(sigma2(j)-h2*k22) - alpha*(sigma3(j)-h2*k23) -...
            theta*(sigma5(j)-h2*k25) - double(j<=(NN+1-taux3) && j>=0) * ...
            (sigma4(j+taux3)-h2*k24)*0.5*(w(j+taux3)+w(j+taux3-1));
        k33=-1+(mu+lamda+delta1+0.5*(v(j)+v(j-1)))*(sigma3(j)-h2*k23) - ...
            lamda*(sigma4(j)-h2*k24) - double(j<=(NN+1-taux2) && j>=0) * ...
            (sigma4(j+taux2)-h2*k24)*0.5*(v(j+taux2)+v(j+taux2-1));
        k34=-1+(mu+gamma+delta2+0.5*(z(j)+z(j-1)))*(sigma4(j)-h2*k24) - ...
            gamma*(sigma5(j)-h2*k25) - double(j<(NN-taux4) && j>0) * ...
            (sigma5(j+taux4)-h2*k25)*0.5*(z(j+taux4)+z(j+taux4-1));
        k35=mu*(sigma5(j)-h2*k25);
        
        
        k41=(mu+beta*E(j-1)/N)*(sigma1(j)-h2*k31)-...
            (beta*(1-u(j-1))*E(j-1)/N)*(sigma2(j)-h2*k32) - beta * ...
            double(j<=(NN+1-taux1) && j>=0) * ...
            (sigma1(j+taux1)-h2*k31)*u(j+taux1-1)*E(j+taux1-1)/N;
        k42=-1+(beta*S(j-1)/N)*(sigma1(j)-h2*k31)-...
            (beta*(1-u(j-1))*S(j-1)/N-mu-alpha-theta-...
            w(j-1))*(sigma2(j)-h2*k32) - alpha*(sigma3(j)-h2*k33) -...
            theta*(sigma5(j)-h2*k35) - double(j<=(NN+1-taux3) && j>=0) * ...
            (sigma4(j+taux3)-h2*k34)*w(j+taux3-1);
        k43=-1+(mu+lamda+delta1+v(j-1))*(sigma3(j)-h2*k33) - ...
            lamda*(sigma4(j)-h2*k34) - double(j<=(NN+1-taux2) && j>=0) * ...
            (sigma4(j+taux2)-h2*k34)*v(j+taux2-1);
        k44=-1+(mu+gamma+delta2+z(j-1))*(sigma4(j)-h2*k34) - ...
            gamma*(sigma5(j)-h2*k35) - double(j<(NN-taux4) && j>0) * ...
            (sigma5(j+taux4)-h2*k35)*z(j+taux4-1);
        k45=mu*(sigma5(j)-h2*k35);
        
        
        sigma1(j-1)=sigma1(j)-(h/6)*(k11+2*k21+2*k31+k41);  % [NN+taux4 x 1]
        sigma2(j-1)=sigma2(j)-(h/6)*(k12+2*k22+2*k32+k42); % [NN+taux4 x 1]
        sigma3(j-1)=sigma3(j)-(h/6)*(k13+2*k23+2*k33+k43); % [NN+taux4 x 1]
        sigma4(j-1)=sigma4(j)-(h/6)*(k14+2*k24+2*k34+k44);  % [NN+taux4 x 1]
        sigma5(j-1)=sigma5(j)-(h/6)*(k15+2*k25+2*k35+k45);  % [NN+taux4 x 1]
        
    end
    
    % sigmas with delays taux
    sigma1_taux1(1:end-taux1) = sigma1(1+taux1 : end); % [(NN+taux4-taux1)x1]
    sigma4_taux2(1:end-taux2) = sigma4(1+taux2 : end); % [(NN+taux4-taux2)x1]
    sigma4_taux3(1:end-taux3)= sigma4(1+taux3 : end); % [(NN+taux4-taux3)x1]
    sigma5_taux4(1:end-taux4) = sigma5(1+taux4 : end); % [(NN+taux4-taux4)x1]
    
    count = count +1;
    
    temp1 = ( sigma2 .* beta .* S .* E - double(tt<=(NN+1-taux1) & tt>=0) .* ...
        sigma1_taux1 .* beta .* S_taux1 .* E_taux1 ) ./(A .* N);
    temp2 = (sigma3 .* I - double(tt<=(NN+1-taux2) & tt>=0) .*...
        sigma4_taux2 .* I_taux2) ./ B;
    temp3 = (sigma2 .* E - double(tt<=(NN+1-taux3) & tt>=0) .* ...
        sigma4_taux3 .* E_taux3 ) ./C;
    temp4 = (sigma4 .* H - double(tt<=(NN+1-taux4) & tt>=0) .*...
        sigma5_taux4 .* H_taux4) ./ D;
    
    u11=min(1,max(0,temp1));
    v1=min(1,max(0,temp2));  w1=min(1,max(0,temp3)); %ancienne config
    z11=min(1,max(0,temp4));
    % v1=0; w1=0; z11=0;
    u=0.5*(u11 + oldu);
    v=0.5*(v1 + oldv); w=0.5*(w1 + oldw);
    z=0.5*(z11 + oldz);
    % v=0; w=0; z=0;
    x=[S E I H R]';
    sigma=[sigma1 sigma2 sigma3 sigma4 sigma5]';
    
    test=delta*sum(abs(u))-sum(abs(oldu-u));
    
    
end
%% code pour enregistrer les valeurs dans fichier covide_retard....mat
tt= tt(1:end-taux4)-1;
z1 = z1(1:end-taux4); z2 = z2(1:end-taux4); z3=z3(1:end-taux4);
z4=z4(1:end-taux4); z5=z5(1:end-taux4);
S = S(1:end-taux4); E = E(1:end-taux4); I=I(1:end-taux4);
H=H(1:end-taux4); R=R(1:end-taux4);
u11 =u11(1:end-taux4); v1 = v1(1:end-taux4); w1 = w1(1:end-taux4);
z11 = z11(1:end-taux4);



end
