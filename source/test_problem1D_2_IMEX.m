clear all
format short
video_on = 0;
%{
Numerical solution to the Gray-Scott ReactionDiffusion equation using Hyperbolic B-spline
non-Linear RD
%}

m = 401
x_l=-20;x_r=20;
len = x_r - x_l;                  
n=2*m;
h=(x_r-x_l)/(m-1);


I2= eye(2);
I= eye(n);
e1=[1 0];
e2=[0 1];

SBP4;
%SBP6_Higher_2018;
%e_1=e_l;e_m=e_r;
%SBP6;
HII =  kron(I2,HI);
L =  [ kron(I2,e_1');kron(I2,e_m')];
%L = sparse([ kron(I2,d1_l);kron(I2,d2_l)]);
P=I-HII*L'*((L*HII*L')\L);
D2 = kron(I2,D2);
D1 = kron(I2,D1);

%D3 = kron(I2,D3);

M = 1/8;
alpha = 1-4*M;
beta = 1+sqrt(alpha) - 2*M;
gamma = (sqrt(2)*(1-3*sqrt(alpha)))/4;
x=linspace(x_l,x_r,m)';	

C1 =  kron([-1 0; 0 1],eye(m)) ;
c2 = ones(n,1)*(-M);
C2 =  diag(c2);
C3 = zeros(n,1); C3(1:m) = M;

V= zeros(n,1) ;          
V(1:m)   = U_exact(x,0,alpha, beta, gamma);
V(m+1:n) =  V_exact(x,0,alpha, beta, gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_1= 50;  
CFL=0.5;% End time
%CFL=0.001;what we found
k=CFL*h
max_itter=floor(t_1/k);
freq = 0;
t = 0;

BC =  zeros(m,1)  ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMEX Runge–Kutta schemes for reaction–diffusion equations
% (2) ESDIRK-ERK 4o6s
embedded_order = 3;
step_reduced = 0.1;          % failed step reduction factor
step_safety = 0.9;          % adaptivity safety factor
step_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);     % coefficients to account for
ONEPSM   = 1+sqrt(eps);     %   floating-point roundoff
Tolerance   = 1.5;             % upper bound on allowed step error
rej = 0;
atol = 1e-7;
rtol = 1e-7;

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
iteration_M = (eye(n)-k*a_s(2,2)*P*D2)\(eye(n));


error_1=0.1 ;error_2=0.1;
K = [0*V 0*V 0*V 0*V 0*V 0*V];
K_hat = [0*V 0*V 0*V 0*V 0*V 0*V];

t1 = cputime;
%for nr_itter=1:10
while t<t_1
    K(:,1) =  D2*V;
    K_hat(:,1) = (C1*uv2(V,m) +C2*V  +C3);
    F_s_total = b_s(1)*K(:,1);  
    G_ns_total = b_ns(1)*K_hat(:,1);
    
    F_s_embedded = b_s_hat(1)*K(:,1);
    G_ns_embedded = b_ns_hat(1)*K_hat(:,1);

    for i=2:6
        BC(1) = U_exact(x_l,t+k*c_s(i), alpha, beta, gamma); BC(m) = U_exact(x_r,t+k*c_s(i), alpha, beta, gamma);
        BC(m+1)= V_exact(x_l,t+k*c_s(i), alpha, beta, gamma); BC(n) = V_exact(x_r,t+k*c_s(i), alpha, beta, gamma);
   
        
        F_s = 0*V;
        G_ns = 0*V;

        for j=1:i-1

          
            F_s = F_s +  a_s(i,j)* K(:,j);
            G_ns = G_ns + a_ns(i,j)* K_hat(:,j);
        end        
           
        W_guess = iteration_M*(P*(V  + k*F_s +k*G_ns)+BC);
   
       
        W_guess = P*(W_guess)+BC;
        K(:,i) = D2*W_guess;
        K_hat(:,i) = (C1*uv2(W_guess,m) +C2*W_guess +C3);

        F_s_total = F_s_total + b_s(i)*K(:,i) ;
        G_ns_total = G_ns_total +b_ns(i)*K_hat(:,i);

        F_s_embedded =  F_s_embedded +b_s_hat(i)*K(:,i);
        G_ns_embedded = G_ns_embedded + b_ns_hat(i)*K_hat(:,i);

       

    end
    %---------------------------------------------------
    
    BC(1) = U_exact(x_l,t+k, alpha, beta, gamma); BC(m) = U_exact(x_r,t+k, alpha, beta, gamma);
    BC(m+1)= V_exact(x_l,t+k, alpha, beta, gamma); BC(n) = V_exact(x_r,t+k, alpha, beta, gamma);
 
 %   V_new = P*(V + k*F_s_total + k*G_ns_total)+BC ;
 %   V = V_new;
 %   t = t+k;
 %Control
    %
    V_new= P*(V + k*F_s_total + k*G_ns_total   )+BC  ;
    V_hat = P*(V+ k*F_s_embedded +k*G_ns_embedded) + BC;
    
    V_error = V_hat-V_new;
    step_error = max(norm(V_error./(rtol*V_new + atol),inf), eps);
    
    %step_error = norm(V_error);%./(V_new ),inf), eps);
    k_old = k;
    if (step_error > Tolerance*ONEPSM*0.1)
        rej = rej + 1;
        k = k*step_reduced;           
    else
        V = V_new; 
        t = t+k;
       %epss = 1e-1;
       epss = 0.1;
       %(PID-controller)
     
       k = step_safety*k*(epss/step_error)^(0.49/embedded_order)*(error_1/epss)^(0.34/embedded_order)*(epss/error_2)^(0.1/embedded_order);

       error_2 = error_1; 
       error_1 = step_error;
       
       %(I-controller)
       %k = step_safety * k * step_error^(-1.0/(embedded_order-1))
    end
   k = min(step_growth*k_old, k) ;
   iteration_M =  (eye(n)-k*a_s(2,2)*P*D2)\(eye(n)) ;

    
    %--------------------------------------------------------
end
t2 = cputime - t1
L_U_2 = sqrt((U_exact(x,t,alpha,beta,gamma)-V(1:m))'*H*(U_exact(x,t,alpha,beta,gamma)-V(1:m)));
L_V_2 = sqrt((V_exact(x,t,alpha,beta,gamma)-V(m+1:2*m))'*H*(V_exact(x,t,alpha,beta,gamma)-V(m+1:2*m)));
err = L_U_2 + L_V_2

%L_U_inf = max(abs(U_exact(x,t,alpha,beta,gamma)-V(1:m)))
%L_V_inf = max(abs(V_exact(x,t,alpha,beta,gamma)-V(m+1:2*m)))


function pq2 = uv2(U,m)
    pq2 = U(1:m).*(U(m+1:2*m)).^2;
    pq2 = [pq2;pq2];
end


function p = U_initial(x, alpha, beta)
    p = (3-sqrt(alpha))/4 - sqrt(2*beta)/4*tanh(sqrt(beta)/4*x);
end

function q = V_initial(x, alpha, beta)
    q = (1+sqrt(alpha))/4 + sqrt(2*beta)/4*tanh(sqrt(beta)/4*x);
end


function p = U_exact(x,t, alpha, beta, gamma)
    p = (3-sqrt(alpha))/4 - sqrt(2*beta)/4*tanh(sqrt(beta)/4*(x-gamma*t));
end

function q = V_exact(x,t, alpha, beta, gamma)
    q = (1+sqrt(alpha))/4 + sqrt(2*beta)/4*tanh(sqrt(beta)/4*(x-gamma*t));
end

function x = ffsolve(varargin)
    [x,~] = fsolve(varargin{:});
end


% New starting point
%https://www.jstor.org/stable/pdf/2158449.pdf?refreqid=excelsior%3A4759cae3f809aadc8e95023d2c00bb04&ab_segments=&origin=&acceptTC=1
