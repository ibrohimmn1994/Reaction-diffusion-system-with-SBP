
%{

non-Linear RD
%}



clear all
format short

method = "RK4" ;
global m


m =101;             %grid pints
x_l=-20;x_r=20;
len = x_r - x_l;                  
n=2*m;
h=(x_r-x_l)/(m-1);


I2=eye(2);
I=eye(n);
e1=[1 0];
e2=[0 1];

SBP6;              % SBP operators
% SBP4

HII = kron(I2,HI);
L = [ kron(I2,e_1');kron(I2,e_m')]; % 
P=I-HII*L'*((L*HII*L')\L);          % projection operator
D2 = kron(I2,D2);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corefficients
M = 1/8;
alpha = 1-4*M;
beta = 1+sqrt(alpha) - 2*M;
gamma = (sqrt(2)*(1-3*sqrt(alpha)))/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(x_l,x_r,m)';	        % grid

%%%%%%%%%%%%%%%%%%%%%%%%%%
C1 = kron([-1 0; 0 1],eye(m)) ;
c2 = ones(n,1)*(-M);
C2 = diag(c2);
C3 = zeros(n,1); C3(1:m) = M;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial data
V= zeros(n,1) ;          
V(1:m)   = U_exact(x,0,alpha, beta, gamma);
V(m+1:n) =  V_exact(x,0,alpha, beta, gamma);
%%%%%%%%%%%%%%%%%%%%%%%%%%
t_1= 2;                    % End time
CFL=0.1;                   
k=CFL*h^2;                 % Time step
max_itter=floor(t_1/k);
t = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%

BC = zeros(m,1)  ;



t_start = cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROCK4
coeff = {D2, C1, C2, C3};

atol =1e-13;
rtol =1e-12;
%[V,t] = ROCK4_solver(0 , t_1 , coeff,0, V, P , BC, k, "analytical_2",x,atol, rtol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMEX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RK4

while t<t_1
    

    W = V ;
    w1 =  (D2*(W)  +C1*uv2(W,m) +C2*W  +C3) ;
    
    BC(1) = U_exact(x_l,t+k/2, alpha, beta, gamma); BC(m) = U_exact(x_r,t+k/2, alpha, beta, gamma);
    BC(m+1)= V_exact(x_l,t+k/2, alpha, beta, gamma); BC(n) = V_exact(x_r,t+k/2, alpha, beta, gamma);


    W = P*(V+k/2*w1) +BC;
    w2 = (D2*(W)  +C1*uv2(W,m) +C2*W  +C3) ;

    BC(1) = U_exact(x_l,t+k/2, alpha, beta, gamma); BC(m) = U_exact(x_r,t+k/2, alpha, beta, gamma);
    BC(m+1)= V_exact(x_l,t+k/2, alpha, beta, gamma); BC(n) = V_exact(x_r,t+k/2, alpha, beta, gamma);

    W = P*(V+k/2*w2) +BC;
    w3 =  (D2*(W)  +C1*uv2(W,m) +C2*W  +C3) ;

    BC(1) = U_exact(x_l,t+k, alpha, beta, gamma); BC(m) = U_exact(x_r,t+k, alpha, beta, gamma);
    BC(m+1)= V_exact(x_l,t+k, alpha, beta, gamma); BC(n) = V_exact(x_r,t+k, alpha, beta, gamma);

    W = P*(V+k*w3) +BC;
    w4 =  (D2*(W)  +C1*uv2(W,m) +C2*W  +C3) ;

    V= P*(V+k/6*(w1+2*w2+2*w3+w4)) +BC;

    t=t+k;
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

t_end = cputime - t_start
L_U_2 = sqrt((U_exact(x,t,alpha,beta,gamma)-V(1:m))'*H*(U_exact(x,t,alpha,beta,gamma)-V(1:m)));
L_V_2 = sqrt((V_exact(x,t,alpha,beta,gamma)-V(m+1:2*m))'*H*(V_exact(x,t,alpha,beta,gamma)-V(m+1:2*m)));

%numerical error
err = L_U_2 + L_V_2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pq2 = uv2(U,m)
    pq2 = U(1:m).*(U(m+1:2*m)).^2;
    pq2 = [pq2;pq2];
end



% analytical solution u-component
function p = U_exact(x,t, alpha, beta, gamma)
    p = (3-sqrt(alpha))/4 - sqrt(2*beta)/4*tanh(sqrt(beta)/4*(x-gamma*t));
end
% analytical solution v-compnent
function q = V_exact(x,t, alpha, beta, gamma)
    q = (1+sqrt(alpha))/4 + sqrt(2*beta)/4*tanh(sqrt(beta)/4*(x-gamma*t));
end





