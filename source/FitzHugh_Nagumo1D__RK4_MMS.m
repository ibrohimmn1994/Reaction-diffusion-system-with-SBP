clear all;
format short


dd1 = 0.005;
dd2 = 0.00125;

u0 = 2;
v0 = 1;

m =501

x_l=0;x_r=1;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

t_1= 1.0;                    % End time
CFL=25;  
k=CFL*h^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


I2=eye(2);
I=eye(m);

SBP6;


%HII = kron(I2,HI);
Lu = [ e_1'; e_m'];
Lv = [ e_1'; e_m'];
Pu=I-HI*Lu'*((Lu*HI*Lu')\Lu);
Pv=I-HI*Lv'*((Lv*HI*Lv')\Lv);


bcu = zeros(2,1);
bcu(1) = u0; bcu(2) = u0;
BCu= HI*Lu'*((Lu*HI*Lu')\bcu);

bcv = zeros(2,1);
bcv(1) = v0; bcv(2) = v0;
BCv= HI*Lv'*((Lv*HI*Lv')\bcv);



max_itter=floor(t_1/k);

x=linspace(x_l,x_r,m)';	
          
%u   = u0 + 8*exp( -((x-0.5)/0.025).^2 ) ;
%v = ones(m,1)*v0;

u = u_exact(x,0);
v = v_exact(x,0);
%coeff = DD*D2;
t=0.0;
freq = 0;

d1 = -0.75;
d2 = 52;
d3 = 10;
d4 = 0.075;
d5 = 0.52;
d6 = 10;
d7 = 9; 
d8 = 180;
u_ =  2.2;
d11 = 0.11;
vt = 46;
d10 = 0.3;
k_ = 0.25 ;
d9 =  140;
tStart = cputime;
i=1;
for nr_itter=1:max_itter

    %_____________ one
    Wu = Pu*u+BCu;
    Wv = Pv*v+BCv;

    [su,sv] = source(x,t);
    [F,G] = reaction(Wu,Wv);

    wu1 =    ( dd1*(D2*Wu)      +F       + su  )  ;
    wv1 =    ( dd2*(D2*Wv)      +G        + sv  ) ;
    %_____________ two

    Wu = Pu*(u+k/2*wu1)+BCu;
    Wv = Pv*(v+k/2*wv1)+BCv;

    [su,sv] = source(x,t+k/2);
    [F,G] = reaction(Wu,Wv);
    
    wu2 =   (dd1*(D2*Wu)         +F     +  su)  ;
    wv2 =   (dd2*(D2*Wv)         +G      + sv )  ;
    %____________ three
    Wu = Pu*(u+k/2*wu2)+BCu;
    Wv = Pv*(v+k/2*wv2)+BCv;

    [su,sv] = source(x,t+k/2);
    [F,G] = reaction(Wu,Wv);

    wu3 =    (dd1*(D2*Wu)        +F    + su ) ;
    wv3 =    (dd2*(D2*Wv)        +G    + sv );
    %_____________ four
    Wu = Pu*(u+k*wu3)+BCu;
    Wv = Pv*(v+k*wv3)+BCv;

    [su,sv] = source(x,t+k);
    [F,G] = reaction(Wu,Wv);

    wu4 =    (dd1*(D2*Wu)        +F    + su  )  ;
    wv4 =    (dd2*(D2*Wv)        +G     + sv  )  ;

    u=u+k/6*(wu1+2*wu2+2*wu3+wu4) ;%+ k* source_u(x,t+k) ;
    v=v+k/6*(wv1+2*wv2+2*wv3+wv4) ;%+ k*source_v(x,t+k) ;
    t=t+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    if t>i*0.1
    t;
    erro_old =err;
    i = i+1
    L_U_2 = sqrt((u_exact(x,t)-u)'*H*(u_exact(x,t)-u));
    L_V_2 = sqrt((v_exact(x,t)-v)'*H*(v_exact(x,t)-v));
    err = L_U_2 + L_V_2
    err/erro_old
    
    end
    %}
    

end

tSCPU = cputime - tStart

L_U_2 = sqrt((u_exact(x,t)-u)'*H*(u_exact(x,t)-u));
L_V_2 = sqrt((v_exact(x,t)-v)'*H*(v_exact(x,t)-v));

err = L_U_2 + L_V_2
disp("#########################")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [su, sv] = source(x,t)

D1 = 0.005;D2 = 0.00125;d1 = -0.75;d2 = 52;d3 = 10;d4 = 0.075;d5 = 0.52;
d6 = 10;d7 = 9; d8 = 180;u_ =  2.2;d11 = 0.11;vt = 46;d10 = 0.3;k = 0.25 ;
d9 =  140;u0 = 2;v0 = 1;

u =  2 + 8*exp( -((x-0.5)/0.025).^2)*exp(t/4);%sin((t+pi)/2);
v =  1 -exp( -((x-0.5)/0.025 ).^2)*(1-1/exp(t));%sin(t/2);
%u = 2 + 8*exp( -((x-0.5)/0.025).^2)*sin((t+pi)/2);
%v = 1 - exp( -((x-0.5)/0.025 ).^2)*sin(t/2);

H_ = heaviside(u-u_);%(Wu-u_)>0;
VV =  58.*log10((u+d7)./d8) ;
hh =  ( 1+tanh( d11.* (VV+vt) ) ).*H_ ;
Vc =   29*  log10(v./(d10-k.*v))   ;
Vk =   58.*log10(u./d9)  ;

G = d4.*hh .* (VV - Vc) + d5.*( 1-exp(-d6.*(v0-v)) );
%differential_part_v = - (cos(t/2)*exp(-(40*x - 20).^2))/2 - D2*(3200*sin(t/2)*exp(-(40*x - 20).^2) - sin(t/2)*exp(-(40*x - 20).^2).*(3200*x - 1600).^2) ;
differential_part_v = -exp(-t)*exp(-(40*x - 20).^2) -D2*(  exp(-(40*x - 20).^2).*(3200*x - 1600).^2*(exp(-t) - 1) - 3200*exp(-(40*x - 20).^2)*(exp(-t) - 1)  );
 

F = d1.*hh .* (VV - Vc) .* (VV - Vk) - d2.*( 1-exp(-d3.*(u-u0)) );
%differential_part_u = +4*cos(t/2 + pi/2)*exp(-(40*x - 20).^2) - D1*(-25600*sin(t/2 + pi/2)*exp(-(40*x - 20).^2) + 8*sin(t/2 + pi/2)*exp(-(40*x - 20).^2).*(3200*x - 1600).^2);
%differential_part_u = -8*exp(-t)*exp(-(40*x - 20).^2) - D1*(8*exp(-t)*exp(-(40*x - 20).^2).*((3200*x - 1600).^2) - 25600*exp(-t)*exp(-(40*x - 20).^2));
%differential_part_u = 4*exp(-(40*x - 20).^2)*exp(t) - D1*(4*exp(-(40*x - 20).^2).*exp(t).*(3200*x - 1600).^2 - 0.5*25600*exp(-(40*x - 20).^2)*exp(t));
differential_part_u = 2*exp(t/4)*exp(-(40*x - 20).^2) - D1*(8*exp(t/4)*exp(-(40*x - 20).^2).*(3200*x - 1600).^2 - 25600*exp(t/4)*exp(-(40*x - 20).^2));
sv = differential_part_v - G; 
su = differential_part_u - F;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, G] = reaction(Wu,Wv)
d1 = -0.75;
d2 = 52;
d3 = 10;
d4 = 0.075;
d5 = 0.52;
d6 = 10;
d7 = 9; 
d8 = 180;
u_ =  2.2;
d11 = 0.11;
vt = 46;
d10 = 0.3;
k_ = 0.25 ;
d9 =  140;
u0 = 2; v0 = 1;

H_ =  heaviside(Wu-u_);%(Wu-u_)>0;
%H_ = (Wu-u_)>0;
VV =  58.*log10((Wu+d7)./d8) ;
hh =  ( 1+tanh( d11.* (VV+vt) ) ).*H_ ;
Vc =   29*  log10(Wv./(d10-k_.*Wv))   ;
Vk =   58.*log10(Wu./d9)  ;
F = d1.*hh .* (VV - Vc) .* (VV - Vk) - d2.*( 1-exp(-d3.*(Wu-u0)) );
G = d4.*hh .* (VV - Vc) + d5.*( 1-exp(-d6.*(v0-Wv)) );
end

function u = u_exact(x,t)

u = 2 + 8*exp( -((x-0.5)/0.025).^2)*exp(t/4);%sin((t+pi)/2)
%u = 2 + 4*exp( -((x-0.5)/0.025).^2)*exp(t)
%u = 2 + 8*exp( -((x-0.5)/0.025).^2)*sin((t+pi)/2);
end

function v = v_exact(x,t)
v = 1 - exp( -((x-0.5)/0.025 ).^2)*(1-1/exp(t));%sin(t/2);
%v = 1 - exp( -((x-0.5)/0.025 ).^2)*sin(t/2);
end



