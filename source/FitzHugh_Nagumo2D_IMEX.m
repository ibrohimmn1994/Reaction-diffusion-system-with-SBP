%{
Finite element analysis of nonlinear reaction–diffusion system of
Fitzhugh–Nagumo type with Robin boundary conditions


%}
clear all
format short
m =101;
x_l=-50.0;x_r=50.0;
len = x_r - x_l;                  
n=m*m;
n_t = 2*n;
h=(x_r-x_l)/(m-1);


CFL=2;
k=CFL*h;
%k = 10^-3
K = 2;
Ik= speye(K);
Im= speye(m);
In= speye(n);
I= speye(K*m*m);
I2 = eye(2);

e1=[1 0];
e2=[0 1];

SBP4;

Dxx = kron( sparse(D2), Im); Dyy = kron(Im, sparse(D2));
HIx = kron( sparse(HI), Im); HIy = kron(Im, sparse(HI));
HI = HIx*HIy;
HI_bar = sparse(kron(Ik, HI)) ;

Hx = kron( sparse(H), Im); Hy = kron(Im,  sparse(H));
H_bar = Hx*Hy;

clear HIx HIy D1 D2 M Q M_U Q_U S_U H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_E = kron( d_m, Im); d_W = kron( d_1, Im);
d_N = kron(Im,  d_m); d_S = kron(Im,  d_1);

Lx = [ kron(I2,d_W); kron(I2,d_E)];
Ly = [kron(I2,d_S) ; kron(I2,d_N)];

Px = I - HI_bar*Lx'*((Lx*HI_bar*Lx')\Lx);
Py = I - HI_bar*Ly'*((Ly*HI_bar*Ly')\Ly);
P = sparse(Px*Py);

Dxx = sparse(kron(Ik,Dxx));
Dyy = sparse(kron(Ik,Dyy));
C1 = ones(n_t,1); C1(n+1:n_t) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
video_on = 0;

t_1= 20;
%t_1 = 2*k;
V = zeros(2*n,1); 
%
V(1:n) =  reshape(u_initial(X,Y),n,1);
V(n+1:2*n) =  zeros(n,1);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0.0;
M = sparse(I-k*a_s(2,2)*P*C1.*(Dxx + Dyy));
iteration_M =  inv(M);%((I-k*a_s(2,2)*P*C1.*(Dxx + Dyy))\I);
tStart = cputime;
while t<t_1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %stage 1
    K1 = C1.*(Dxx*V + Dyy*V);
    K_hat1 = (R(V,n) + S(X,Y,t,n));

    F_s_total = b_s(1)*K1;  
    G_ns_total = b_ns(1)*K_hat1;
    
  %  F_s_embedded = b_s_hat(1)*K1;
  %  G_ns_embedded = b_ns_hat(1)*K_hat1;
    %-----------------------------------------------
    %stage 2
    F_s =  a_s(2,1)*K1;
    G_ns = a_ns(2,1)*K_hat1;

    W_guess = P*( iteration_M*(P*(V  + k*F_s + k*G_ns)) );

    K2 = C1.*(Dxx*W_guess + Dyy*W_guess);
    K_hat2 = (R(W_guess,n) + S(X,Y,t+k*c_s(2),n));

    F_s_total = F_s_total + b_s(2)*K2;  
    G_ns_total = G_ns_total + b_ns(2)*K_hat2;
    
 %   F_s_embedded = F_s_embedded +b_s_hat(2)*K2;
 %   G_ns_embedded = G_ns_embedded + b_ns_hat(2)*K_hat2;
    %--------------------------------------------------------
    % stage 3

    F_s =  a_s(3,1)*K1 + a_s(3,2)*K2;
    G_ns = a_ns(3,1)*K_hat1 + a_ns(3,2)*K_hat2;

    W_guess = P*( iteration_M*(P*(V  + k*F_s + k*G_ns)) );

    K3 = C1.*(Dxx*W_guess + Dyy*W_guess);
    K_hat3 = (R(W_guess,n) + S(X,Y,t+k*c_s(3),n));

    F_s_total = F_s_total + b_s(3)*K3;  
    G_ns_total = G_ns_total + b_ns(3)*K_hat3;
    
 %   F_s_embedded = F_s_embedded +b_s_hat(3)*K3;
 %   G_ns_embedded = G_ns_embedded + b_ns_hat(3)*K_hat3;
    
    %--------------------------------------------------------
    % stage 4
    F_s =  a_s(4,1)*K1 + a_s(4,2)*K2 + a_s(4,3)*K3;
    G_ns = a_ns(4,1)*K_hat1 + a_ns(4,2)*K_hat2 + a_ns(4,3)*K_hat3;

    W_guess = P*( iteration_M*(P*(V  + k*F_s + k*G_ns)) );

    K4 = C1.*(Dxx*W_guess + Dyy*W_guess);
    K_hat4 = (R(W_guess,n) + S(X,Y,t+k*c_s(4),n));

    F_s_total = F_s_total + b_s(4)*K4;  
    G_ns_total = G_ns_total + b_ns(4)*K_hat4;
    
 %   F_s_embedded = F_s_embedded +b_s_hat(4)*K4;
 %   G_ns_embedded = G_ns_embedded + b_ns_hat(4)*K_hat4;
    
    %--------------------------------------------------------
    % stage 5

    F_s =  a_s(5,1)*K1 + a_s(5,2)*K2 + a_s(5,3)*K3 +a_s(5,4)*K4;
    G_ns = a_ns(5,1)*K_hat1 + a_ns(5,2)*K_hat2 + a_ns(5,3)*K_hat3 + a_ns(5,4)*K_hat4;

    W_guess = P*( iteration_M*(P*(V  + k*F_s + k*G_ns)) );

    K5 = C1.*(Dxx*W_guess + Dyy*W_guess);
    K_hat5 = (R(W_guess,n) + S(X,Y,t+k*c_s(5),n));

    F_s_total = F_s_total + b_s(5)*K5;  
    G_ns_total = G_ns_total + b_ns(5)*K_hat5;
    
  %  F_s_embedded = F_s_embedded +b_s_hat(5)*K5;
  %  G_ns_embedded = G_ns_embedded + b_ns_hat(5)*K_hat5;
    

    %--------------------------------------------------------
    % stage 6
    F_s =  a_s(6,1)*K1 + a_s(6,2)*K2 + a_s(6,3)*K3 +a_s(6,4)*K4 + a_s(6,5)*K5;
    G_ns = a_ns(6,1)*K_hat1 + a_ns(6,2)*K_hat2 + a_ns(6,3)*K_hat3 + a_ns(6,4)*K_hat4 + a_ns(6,5)*K_hat5;

    W_guess = P*( iteration_M*(P*(V  + k*F_s + k*G_ns)) );

    K6 = C1.*(Dxx*W_guess + Dyy*W_guess);
    K_hat6 = (R(W_guess,n) + S(X,Y,t+k*c_s(6),n));

    F_s_total = F_s_total + b_s(6)*K6;  
    G_ns_total = G_ns_total + b_ns(6)*K_hat6;
    
  %  F_s_embedded = F_s_embedded +b_s_hat(6)*K6;
  %  G_ns_embedded = G_ns_embedded + b_ns_hat(6)*K_hat6;

    V= P*(V + k*F_s_total + k*G_ns_total   )  ;
   %
    t = t+k;
    %--------------------------------------------------------
    %Control
    %{
    V_new= P*(V + k*F_s_total + k*G_ns_total   )  ;
    V_hat = P*(V+ k*F_s_embedded +k*G_ns_embedded) ;
    
    V_error = V_hat-V_new;
    step_error = max(norm(V_error./(rtol*V_new + atol),inf), eps);
    
    k_old = k;
    if (step_error > Tolerance*ONEPSM*0.1)
        rej = rej + 1;
        k = k*step_reduced;           
    else
        V = V_new; 
        t = t+k
       
     
       %(PID-controller)
       epss = 0.1;
       k = step_safety*k*(epss/step_error)^(0.49/embedded_order)*(error_1/epss)^(0.34/embedded_order)*(epss/error_2)^(0.1/embedded_order);
       error_2 = error_1; 
       error_1 = step_error;
       
       %(I-controller)
       %k = step_safety * k * step_error^(-1.0/(embedded_order-1))
    end
   k = min(step_growth*k_old, k) ;

   iteration_M =  ((I-k*a_s(2,2)*P*C1.*(Dxx + Dyy))\I);
    %}

end
tEnd = cputime - tStart
%

%
U_approx = reshape(V(1:n),n,1);
V_approx = reshape(V(n+1:end),n,1) ;
U_exact = reshape(u_exact(X,Y,t),n,1);
V_exact = reshape(v_exact(X,Y,t),n,1);

L_U_2 = sqrt((U_exact-U_approx)'*H_bar*(U_exact-U_approx));
L_V_2 = sqrt((V_exact-V_approx)'*H_bar*(V_exact-V_approx));
err = L_U_2 + L_V_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function ui = u_initial(x,y)
ui = exp(- ((x).^2 + (y).^2)./16);
end

function u = u_exact(x,y,t)
u  = exp(- ((x).^2 + (y).^2)./16)*(1-tanh(t/100));
end

function v = v_exact(x,y,t)
v  = exp(- ((x).^2 + (y).^2)./16)*(tanh(t/100));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = S(x,y,t,n)
u  = exp(- ((x).^2 + (y).^2)./16)*(1-tanh(t/100));
v  = exp(- ((x).^2 + (y).^2)./16)*(tanh(t/100));

dudt = exp(- x.^2/16 - y.^2/16)*(tanh(t/100)^2/100 - 1/100);
dudx2 = (exp(- x.^2/16 - y.^2/16)*(tanh(t/100) - 1))/8 - ((x.^2).*exp(- x.^2/16 - y.^2/16)*(tanh(t/100) - 1))/64;
dudy2 = (exp(- x.^2/16 - y.^2/16)*(tanh(t/100) - 1))/8 - ((y.^2).*exp(- x.^2/16 - y.^2/16)*(tanh(t/100) - 1))/64;
dvdt = -exp(- x.^2/16 - y.^2/16)*(tanh(t/100)^2/100 - 1/100);

dudt = reshape(dudt,n,1);
dudx2 = reshape(dudx2,n,1);
dudy2 = reshape(dudy2,n,1);
dvdt = reshape(dvdt,n,1);
u = reshape(u,n,1);
v = reshape(v,n,1);
su = dudt - dudx2 - dudy2 - (-u.*(u-1).*(u-0.15) - v);
sv = dvdt - (0.005*u - 0.025*v);
s = [su; sv];
end
%
function r = R(W,n)
u = W(1:n);
v = W(n+1:2*n);
f = -u.*(u-1).*(u-0.15) - v;
g = 0.005*u - 0.025*v;
r = [f;g];
end
%{
function X = CG_inv(A)
tol = 1e-3;
n=length(A);
X = speye(n);
B = speye(n);
R = B - A*X;
P = R; 
for k = 1:100000
    W = A*P;
    alpha =sum(sum(R'*R))/sum(sum(P'*W));
    X = X + alpha*P;      
    R_old=R;
    R = R - alpha*W;      
    if normest(R)<tol%(sqrt(sum(abs(R(:)).^2))<tol)
        %(sqrt(sum(abs(R(:)).^2)))<tol%( norm(R) < tol )
       break;
     end
       
     beta = sum(sum(R'*R))/sum(sum(R_old'*R_old));
     P = R + beta*P; 
end  

end
%}

function X = CG0(A,b)
    tol = 1e-3;
    n = length(A);
    x = sparse(zeros(n,1));    
    r = b - A*x;
    p = r; 
    for k = 1:n      
      w = A*p;
      alpha =(r'*r)/(p'*w);
      x = x + alpha*p;      
      r_old=r;
      r = r - alpha*w;      
      if( normest(r) < tol )
        break;
      end     
      B = (r'*r)/(r_old'*r_old);
      p = r + B*p; 
    end
    X=x;
end


