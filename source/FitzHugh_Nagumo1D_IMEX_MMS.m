clear all;
format short


dd1 = 0.005;
dd2 = 0.00125;

u0 = 2;
v0 = 1;

m =301;

x_l=0;x_r=1;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

t_1= 1;                           % End time
CFL = 0.5; 
k=CFL*h^1;

%########################################################################
%IMEX Runge–Kutta schemes for reaction–diffusion equations
% (2) ESDIRK-ERK 4o6s

a_ns = [0 0 0 0 0 0;
    1/2 0 0 0 0 0;
    13861/62500 6889/62500 0 0 0 0;
    -116923316275/2393684061468 -2731218467317/15368042101831 9408046702089/11113171139209 0 0 0;
    -451086348788/2902428689909 -2682348792572/7519795681897 12662868775082/11960479115383 3355817975965/11060851509271 0 0;
    647845179188/3216320057751 73281519250/8382639484533 552539513391/3454668386233 3354512671639/8306763924573 4040/17871 0];

b_ns = [82889/524892 0 15625/83664 69875/102672 -2260/8211 1/4];
b_ns_hat = [4586570599/29645900160 0 178811875/945068544 814220225/1159782912 -3700637/11593932 61727/225920];

a_s = [0 0 0 0 0 0 ;
    1/4 1/4 0 0 0 0 ;
    8611/62500 -1743/31250 1/4 0 0 0 ;
    5012029/34652500 -654441/2922500 174375/388108 1/4 0 0;
    15267082809/155376265600 -71443401/120774400 730878875/902184768 2285395/8070912 1/4 0;
    82889/524892 0 15625/83664 69875/102672 -2260/8211 1/4];

b_s = [82889/524892 0 15625/83664 69875/102672 -2260/8211 1/4];
b_s_hat =[4586570599/29645900160 0 178811875/945068544 814220225/1159782912 -3700637/11593932 61727/225920];

c_s = [0 1/2 83/250 31/50 17/20 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I2=eye(2);
I=eye(2*m);

SBP4;


HII = kron(I2,HI);
L = [ kron(I2,e_1'); kron(I2,e_m')];
P=I-HII*L'*((L*HII*L')\L);
bc = zeros(4,1);
bc(1) = u0; bc(3) = u0;
bc(2)=v0; bc(4) = v0;
BC = HII*L'*((L*HII*L')\bc);

D2 = kron(I2,D2);

max_itter=floor(t_1/k);

x=linspace(x_l,x_r,m)';	
          
u   = u0 + 8*exp( -((x-0.5)/0.025).^2 ) ;
v = ones(m,1)*v0;
V = [u;v];

dd = ones(n,1)';
dd(1:m) = dd1 ;
dd(m+1:n) = dd2;
DD =   diag(dd)  ;
%coeff = DD*D2;
t=0.0;
freq = 0;


iteration_M =  (eye(n)-k*a_s(2,2)*P*DD*D2)\(eye(n)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W = [ 0*V 0*V 0*V 0*V 0*V 0*V];
K = [0*V 0*V 0*V 0*V 0*V 0*V];
K_hat = [0*V 0*V 0*V 0*V 0*V 0*V];
tStart = cputime;
for nr_itter=1:max_itter

   % W(:,1) = V;
    K(:,1) = DD*D2*V;
    K_hat(:,1) = (R(V,m) + S(x,t));
    F_s_total = b_s(1)*K(:,1);  
    G_ns_total = b_ns(1)*K_hat(:,1);
    
    

    for i=2:6
        F_s = 0*V;
        G_ns = 0*V;

        for j=1:i-1    
            F_s = F_s + a_s(i,j) * K(:,j);
            G_ns = G_ns +a_ns(i,j) * K_hat(:,j);
        end        
           
        W_guess = iteration_M*(P*(V  + k*F_s + k*G_ns)+BC);
       
        W_guess = P*(W_guess)+BC;

        K(:,i) = DD*D2*W_guess;
        K_hat(:,i) = R(W_guess,m) +S(x,t+k*c_s(i));

        F_s_total = F_s_total + b_s(i)*K(:,i) ;
        G_ns_total = G_ns_total +b_ns(i)*K_hat(:,i);

    end
    %--------------------------------------------------
    V = P*(V + k*F_s_total + k*G_ns_total   )+BC  ;
    t = t+k;

end
tCPU = cputime - tStart
V = P*V+BC;

u_aprox = V(1:m);
v_aprox = V(m+1:2*m);

L_U_2 = sqrt((u_exact(x,t)-u_aprox)'*H*(u_exact(x,t)-u_aprox));
L_V_2 = sqrt((v_exact(x,t)-v_aprox)'*H*(v_exact(x,t)-v_aprox));


err = L_U_2 + L_V_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Source terms that we are gonna work with

function reaction = R(W,m)

d1 = -0.75;d2 = 52;d3 = 10;d4 = 0.075;d5 = 0.52;
d6 = 10;d7 = 9; d8 = 180;u_ =  2.2;d11 = 0.11;vt = 46;d10 = 0.3;k_ = 0.25 ;
d9 =  140;u0 = 2;v0 = 1;

Wu = W(1:m);
Wv= W(m+1:2*m);
H_ = heaviside(Wu-u_);%(Wu-u_)>0;
VV =  58.*log10((Wu+d7)./d8) ;
hh =  ( 1+tanh( d11.* (VV+vt) ) ).*H_ ;
Vc =   29*  log10(Wv./(d10-k_.*Wv))   ;
Vk =   58.*log10(Wu./d9)  ;
F = d1.*hh .* (VV - Vc) .* (VV - Vk) - d2.*( 1-exp(-d3.*(Wu-u0)) );
G = d4.*hh .* (VV - Vc) + d5.*( 1-exp(-d6.*(v0-Wv)) );

reaction = [F;G];
end


function source = S(x,t)
D1 = 0.005;D2 = 0.00125;d1 = -0.75;d2 = 52;d3 = 10;d4 = 0.075;d5 = 0.52;
d6 = 10;d7 = 9; d8 = 180;u_ =  2.2;d11 = 0.11;vt = 46;d10 = 0.3;k = 0.25 ;
d9 =  140;u0 = 2;v0 = 1;

%u = 2 + 8*exp( -((x-0.5)/0.025).^2)*sin((t+pi)/2);
%v = 1 - exp( -((x-0.5)/0.025 ).^2)*sin(t/2);
u =  2 + 8*exp( -((x-0.5)/0.025).^2)*exp(t/4);
v =  1 -exp( -((x-0.5)/0.025 ).^2)*(1-1/exp(t));

H_ = heaviside(u-u_);%(Wu-u_)>0;
VV =  58.*log10((u+d7)./d8) ;
hh =  ( 1+tanh( d11.* (VV+vt) ) ).*H_ ;
Vc =   29*  log10(v./(d10-k.*v))   ;
Vk =   58.*log10(u./d9)  ;


F = d1.*hh .* (VV - Vc) .* (VV - Vk) - d2.*( 1-exp(-d3.*(u-u0)) );
%differential_part_u = 4*cos(t/2 + pi/2)*exp(-(40*x - 20).^2) + D1*(25600*sin(t/2 + pi/2)*exp(-(40*x - 20).^2) - 8*sin(t/2 + pi/2)*exp(-(40*x - 20).^2).*(3200*x - 1600).^2);
differential_part_u = 2*exp(t/4)*exp(-(40*x - 20).^2) - D1*(8*exp(t/4)*exp(-(40*x - 20).^2).*(3200*x - 1600).^2 - 25600*exp(t/4)*exp(-(40*x - 20).^2));

su = differential_part_u - F;

G = d4.*hh .* (VV - Vc) + d5.*( 1-exp(-d6.*(v0-v)) );
%differential_part_v = - D2*(3200*sin(t/2)*exp(-(40*x - 20).^2) - sin(t/2)*exp(-(40*x - 20).^2).*(3200*x - 1600).^2) - (cos(t/2)*exp(-(40*x - 20).^2))/2;
differential_part_v = -exp(-t)*exp(-(40*x - 20).^2) -D2*(  exp(-(40*x - 20).^2).*(3200*x - 1600).^2*(exp(-t) - 1) - 3200*exp(-(40*x - 20).^2)*(exp(-t) - 1)  );

sv = differential_part_v - G;

source = [su;sv];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = u_exact(x,t)
u =  2 + 8*exp( -((x-0.5)/0.025).^2)*exp(t/4);

end

function v = v_exact(x,t)
v =  1 -exp( -((x-0.5)/0.025 ).^2)*(1-1/exp(t));
end
