%{
Finite element analysis of nonlinear reaction–diffusion system of
Fitzhugh–Nagumo type with Robin boundary conditions


%}

atol = 1e-8;
rtol = 1e-7;
m = 401
x_l=-50.0;x_r=50.0;
len = x_r - x_l;                  
n=m*m;
n_t = 2*n;
h=(x_r-x_l)/(m-1);


CFL=0.5;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear HIx HIy D1 D2 M Q M_U Q_U S_U H;
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

V = zeros(2*n,1); 
%
V(1:n) =  reshape(u_initial(X,Y),n,1);
V(n+1:2*n) =  zeros(n,1);
V = sparse(V);
%
%__________________________________________________________________________

coeff = {C1, Dxx, Dyy};
x = {X,Y};
BC = 0;


tStart = cputime;
[V,t] = ROCK4_solver(0 , t_1 , coeff, V, P , BC, k, "FHN",x,atol,rtol);

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0.0;
while t<t_1
    % stage1
    W = V ;
    
    w1 =   C1.*(Dxx*W + Dyy*W) +  R(W,n) + S(X,Y,t,n);% + C2.*W + C3*W;   %(D2*(W)  +C1*uv2(W,m) +C2*W  +C3) ;
    %______________________________________________________________________
    %stage2
    

    W = P*(V+k/2*w1) ;
    w2 = C1.*(Dxx*W + Dyy*W) + R(W,n) +S(X,Y,t+k/2,n);%+ C2.*W + C3*W;   
    %______________________________________________________________________
    %stage3
    
   
    W = P*(V+k/2*w2) ;
    w3 =  C1.*(Dxx*W + Dyy*W) + R(W,n) + S(X,Y,t+k/2,n);% + C2.*W + C3*W;   
   
    %______________________________________________________________________
    %stage 4
    

    W = P*(V+k*w3) ;
    w4 =  C1.*(Dxx*W + Dyy*W) + R(W,n) + S(X,Y,t+k,n);%+ C2.*W+ C3*W;   
    %_____________________________________________________________________
    V= P*(V+k/6*(w1+2*w2+2*w3+w4)) ;
    
    
    t=t+k;

    %{
    if (mod(freq,1000)==0)
        t
        U_approx = reshape(V(1:n),n,1);
        V_approx = reshape(V(n+1:end),n,1) ;
        U_exact = reshape(u_exact(X,Y,t),n,1);
        V_exact = reshape(v_exact(X,Y,t),n,1);

        L_U_2 = sqrt((U_exact-U_approx)'*H_bar*(U_exact-U_approx));
        L_V_2 = sqrt((V_exact-V_approx)'*H_bar*(V_exact-V_approx));
        err = L_U_2 + L_V_2
    end
    %}
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tEnd = cputime - tStart
%

U_approx = reshape(V(1:n),n,1);
V_approx = reshape(V(n+1:end),n,1) ;
U_exact = reshape(u_exact(X,Y,t),n,1);
V_exact = reshape(v_exact(X,Y,t),n,1);

L_U_2 = sqrt((U_exact-U_approx)'*H_bar*(U_exact-U_approx));
L_V_2 = sqrt((V_exact-V_approx)'*H_bar*(V_exact-V_approx));
err = L_U_2 + L_V_2


%}
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

