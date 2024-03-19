clear all;
format short
global m
global u0
global v0
DD1 = 0.005;
DD2 = 0.00125;

u0 = 2;
v0 = 1;



video_on = 0;



m =151 ;
atol = 1e-7;%1e-4;
rtol = 1e-7;%1e-3;

x_l=0;x_r=1;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);

t_l= 1.0;%0.64;%1.0768;                    % End time
CFL=0.5;
k= CFL*h;
%k =0.200*h;
k;



I2=eye(2);
I=eye(n);

SBP6;


HII = kron(I2,HI);
L = [ kron(I2,e_1'); kron(I2,e_m')];
P=I-HII*L'*((L*HII*L')\L);

bc = zeros(4,1);
bc(1) = u0; bc(3) = u0;
bc(2)=v0; bc(4) = v0;
BC = HII*L'*((L*HII*L')\bc);

          
D2 = kron(I2,D2);

max_itter=floor(t_l/k);

x=linspace(x_l,x_r,m)';	
V=zeros(n,1) ;          
V(1:m)   = u0 + 8*exp( -((x-0.5)/0.025).^2 );%+ 8*exp( -((x-2.5)/0.025).^2 );
V(m+1:n) = ones(m,1)*v0;

dd = ones(n,1)';
dd(1:m) = DD1 ;
dd(m+1:n) = DD2;
DD =  diag(dd)  ;
t=0.0;
freq = 0;
coeff = DD*D2;

%{
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
%}
tStart = cputime;
%tic;
[V,t_l] = ROCK4_solver(0 , t_l , coeff, V, P , BC, k, "solitons",x,m,H,atol,rtol);


tCPU = cputime - tStart
%toc


u1 = V(1:m);
v1 = V(m+1:2*m);
L_U_2 = sqrt((FHN_1D_MMS.u_exact(x,t_l)-u1)'*H*(FHN_1D_MMS.u_exact(x,t_l)-u1));
L_V_2 = sqrt((FHN_1D_MMS.v_exact(x,t_l)-v1)'*H*(FHN_1D_MMS.v_exact(x,t_l)-v1));
err = L_U_2 + L_V_2

%plot(x,V(1:m),'r','LineWidth',1);


