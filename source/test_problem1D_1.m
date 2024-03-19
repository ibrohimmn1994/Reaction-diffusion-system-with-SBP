clear all
%{

linear RD
%}
%------------------------------------
global m
m = 501
x_l=0;x_r=pi/2;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

I2 = eye(2);
I =  eye(n);
e1=[1 0];
e2=[0 1];

SBP6;%<<-------

HII =  kron(I2,HI);
L =  [ kron(I2,d_1) ; kron(I2,e_m') ];
%L = sparse([ kron(e1,d_1) ; kron(e1,e_m');kron(e2,d_1) ; kron(e2,e_m') ]);
P=I-HII*L'*((L*HII*L')\L);
D2 = kron(I2,D2);

%------------------------------------
% We don't really need this as the boundary has zero values, it will be
% just adding zeros
%{
bc = zeros(4,1);
bc(1) = 0; bc(2) = 0;
bc(3) = 0; bc(4) = 0;
BC=  (HII*L'*(L*HII*L'))\bc;
%}
%------------------------------------


x=linspace(x_l,x_r,m)';	
t_1= 1.0;                    % End time
CFL=1;
k = CFL*h^2;
max_itter=floor(t_1/k);
t = 0;

%
%reaction dominant
%{
a = 0.1;
b = 0.01;
d = 1;
%}
%diffusion dominant
%
a = 2;
b = 1;
d = 0.001;
%

V=sparse(zeros(n,1) );          
V(1:m)   = U_exact(x,0,a,b,d);
V(m+1:n) =  V_exact(x,0,a,b,d);

C1 = sparse(  diag(ones(n,1)*d)  );
c2 = zeros(n,1);
c2(1:m) = -a; c2(m+1:2*m) = -b;
C2 = sparse(diag(c2));

C3 = sparse(kron([0 1;0 0],eye(m)) );

coeff = {C1*D2, C2, C3};
%P = 1;
%############################################################################
% Related to IMEX
embedded_order = 3;
step_reduced = 0.1;          % failed step reduction factor
step_safety = 0.9;          % adaptivity safety factor
step_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);     % coefficients to account for
ONEPSM   = 1+sqrt(eps);     %   floating-point roundoff
Tolerance   = 1.5;             % upper bound on allowed step error
rej = 0;
atol = 1e-9;
rtol = 1e-9;

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
t1 = cputime;
%#############################################################################
%RK4 part
for nr_itter=1:max_itter
   
    W = P*V ;
    w1 =  P*(C1*D2*(W)  +C2*W +C3*W ) ;

    W = P*(V+k/2*w1) ;
    w2 =  P*(C1*D2*(W)  +C2*W +C3*W ) ;

    W = P*(V+k/2*w2) ;
    w3 =  P*(C1*D2*(W)  +C2*W +C3*W ) ;

    W = P*(V+k*w3) ;
    w4 =  P*(C1*D2*(W)  +C2*W +C3*W ) ;
  
    V=V+k/6*(w1+2*w2+2*w3+w4);
  
    t=t+k;   
end
%IMEX part
%############################################################################
%{
error_0 = 1 ;error_1=1 ;error_2=1;
iteration_M =  (eye(n)-k*a_s(2,2)*P*C1*D2)\(eye(n)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = [0*V 0*V 0*V 0*V 0*V 0*V];
K_hat = [0*V 0*V 0*V 0*V 0*V 0*V];
tStart = cputime;
for nr_itter=1:max_itter

   % W(:,1) = V;
    K(:,1) = C1*D2*V;
    K_hat(:,1) = C2*V +C3*V ;
    F_s_total = b_s(1)*K(:,1);  
    G_ns_total = b_ns(1)*K_hat(:,1);
    
    F_s_embedded = b_s_hat(1)*K(:,1);
    G_ns_embedded = b_ns_hat(1)*K_hat(:,1);

    for i=2:6
        F_s = 0*V;
        G_ns = 0*V;

        for j=1:i-1    
            F_s = F_s + a_s(i,j) * K(:,j);
            G_ns = G_ns +a_ns(i,j) * K_hat(:,j);
        end        
           
        W_guess = iteration_M*(P*(V  + k*F_s + k*G_ns));
       
        W_guess = P*(W_guess);

        K(:,i) = C1*D2*W_guess;
        K_hat(:,i) = C2*W_guess +C3*W_guess ;

        F_s_total = F_s_total + b_s(i)*K(:,i) ;
        G_ns_total = G_ns_total +b_ns(i)*K_hat(:,i);

    end
    %--------------------------------------------------
    V_new = P*(V + k*F_s_total + k*G_ns_total   )  ;
   % t = t+k;
    %--------------------------------------------------------
    %Control
    V_hat = P*(V+ k*F_s_embedded +k*G_ns_embedded) ;
    V_error = V_hat-V_new;
    step_error = max(norm(V_error./(rtol*V_new + atol),inf), eps);
    k_old = k;
    if (step_error > Tolerance*ONEPSM*0.1)
        rej = rej + 1;
        k = k*step_reduced;           
    else
        V = V_new;                      
       %epss = 1e-1;
       epss = 0.001;
       %(PID-controller)
       k = step_safety*k*(epss/step_error)^(0.49/embedded_order)*(error_1/epss)^(0.34/embedded_order)*(epss/error_2)^(0.1/embedded_order);
       error_2 = error_1; 
       error_1 = step_error;
       %(I-controller)
       %k = step_safety * k * step_error^(-1.0/(embedded_order));
    end
    k = min(step_growth*k_old, k);   
   

end
%}       
%###########################################################################
t2 = cputime - t1

L_U_2 = sqrt((U_exact(x,t,a,b,d)-V(1:m))'*H*(U_exact(x,t,a,b,d)-V(1:m)));
L_V_2 = sqrt((V_exact(x,t,a,b,d)-V(m+1:2*m))'*H*(V_exact(x,t,a,b,d)-V(m+1:2*m)));
err = L_U_2 + L_V_2

%############################################################################
function u = U_exact(x,t,a,b,d)
    u = ( exp(-(a+d)*t) + exp(-(b+d)*t) )*cos(x);
end
function v = V_exact(x,t,a,b,d)   
    v = (a-b)*( exp(-(b+d)*t) )*cos(x);
end
%##############################################################################