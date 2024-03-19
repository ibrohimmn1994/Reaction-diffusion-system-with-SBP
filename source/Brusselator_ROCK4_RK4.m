


atol=1e-8;
rtol=1e-8;
%Time_integrator = "ROCK4";
Time_integrator = "RK4";


t_1= 40;
m = 25%3.9150e-04
x_l=0.0;x_r=1.0;
len = x_r - x_l;                  
n=m*m;
n_t = 2*n;
h=(x_r-x_l)/(m-1);


%CFL=1;
CFL =0.01;
k=CFL*h^2;

K = 2;
Ik=speye(K);
Im=speye(m);
In=speye(n);
I=speye(K*m*m);
I2 = speye(2);
e1=[1 0];
e2=[0 1];

SBP4;



Dxx = kron(sparse(D2), Im); Dyy = kron(Im, sparse(D2));
HIx = kron(sparse(HI), Im); HIy = kron(Im, sparse(HI));
HI = HIx*HIy;
HI_bar = sparse(kron(Ik, HI)) ;

Hx = kron(sparse(H), Im); Hy = kron(Im, sparse(H));
H_bar = Hx*Hy;
clear HIx HIy D1 M Q M_U Q_U S_U H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_E = kron(sparse(e_m), Im); e_W = kron(sparse(e_1), Im);
e_N = kron(Im, sparse(e_m)); e_S = kron(Im, sparse(e_1));

Lx = [kron(I2,e_W');  kron(I2,e_E') ];
Ly = [kron(I2,e_S') ; kron(I2,e_N') ];


Px = I - HI_bar*Lx'*((Lx*HI_bar*Lx')\Lx);
Py = I - HI_bar*Ly'*((Ly*HI_bar*Ly')\Ly);
P = Px*Py;
bcx = zeros(4*m,1);bcy = zeros(4*m,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


boundary_matrix_x = sparse( HI_bar*Lx'*((Lx*HI_bar*Lx')\eye(4*m)) );
boundary_matrix_y = sparse( HI_bar*Ly'*((Ly*HI_bar*Ly')\eye(4*m)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear e1 e2 e3 e_m e_1;

CC = ones(n,1);
CC(1)=0.5; CC(m)=0.5; CC(end)=0.5; CC(end-m+1)=0.5;
CC = [CC;CC];
%BC = (BCx+BCy).*CC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxx = kron(Ik,Dxx);
Dyy = kron(Ik,Dyy);

C1 = 0.25;
C2 = ones(n_t,1); C2(n+1:n_t) = -1;
C3 = ones(n_t,1); C3(1:n) = -2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


V = zeros(2*n,1);        
V(1:n) =  reshape(u_exact(X,Y,0),n,1);
V(n+1:2*n) =  reshape(v_exact(X,Y,0),n,1);
V = sparse(V);
%__________________________________________________________________________
% ROCK4 stuff
BC_coeff = {boundary_matrix_x, boundary_matrix_y, CC};
DE_coeff = {C1, Dxx, Dyy, C2, C3};

BC = {bcx, bcy};

if Time_integrator == "ROCK4"
tStart = cputime;

[V, t] = ROCK4_solver(0, t_1,DE_coeff,BC_coeff, V, P , BC, k, "Brusselator",x,atol,rtol);
tEnd = cputime - tStart
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Time_integrator == "RK4"
tStart = cputime;
%
t=0.0;
while t< t_1
    % stage1
    W = V ;
  
    w1 =   C1*(Dxx*W + Dyy*W) + C2.*U2V(W,n) + C3.*[W(1:n);W(1:n)];   %(D2*(W)  +C1*uv2(W,m) +C2*W  +C3) ;
    %______________________________________________________________________
    %stage2
    [bcx, bcy] = boundary(bcx,bcy,x,t+k/2);
   
    BCx = boundary_matrix_x*bcx;

   
    BCy = boundary_matrix_y*bcy;
    
    BC = (BCx+BCy).*CC;

    W = P*(V+k/2*w1) +BC;
    w2 = C1*(Dxx*W + Dyy*W) + C2.*U2V(W,n) + C3.*[W(1:n);W(1:n)];   
    %______________________________________________________________________
    %stage3
    [bcx, bcy] = boundary(bcx,bcy,x,t+k/2);

    BCx = boundary_matrix_x*bcx;

   
    BCy = boundary_matrix_y*bcy;
    
    BC = (BCx+BCy).*CC;
   
    W = P*(V+k/2*w2) +BC;
    w3 =  C1*(Dxx*W + Dyy*W) + C2.*U2V(W,n) + C3.*[W(1:n);W(1:n)];   
   
    %______________________________________________________________________
    %stage 4
    [bcx, bcy] = boundary(bcx,bcy,x,t+k);
   
    BCx = boundary_matrix_x*bcx;

 
    BCy = boundary_matrix_y*bcy;
    
    BC = (BCx+BCy).*CC;

    W = P*(V+k*w3) +BC;
    w4 =  C1*(Dxx*W + Dyy*W) + C2.*U2V(W,n) + C3.*[W(1:n);W(1:n)];   
    %_____________________________________________________________________
    V= P*(V+k/6*(w1+2*w2+2*w3+w4)) +BC;
    
    
    t=t+k;
end
%
tEnd = cputime - tStart
end

U_approx = reshape(V(1:n),n,1);
V_approx = reshape(V(n+1:end),n,1) ;
U_exact = reshape(u_exact(X,Y,t),n,1);
V_exact = reshape(v_exact(X,Y,t),n,1);

L_U_2 = sqrt((U_exact-U_approx)'*H_bar*(U_exact-U_approx));
L_V_2 = sqrt((V_exact-V_approx)'*H_bar*(V_exact-V_approx));
err = L_U_2 + L_V_2




%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ui = u_exact(x,y,t)
ui = exp(-x-y-0.5*t);
end
function vi = v_exact(x,y,t)
vi = exp(+x+y+0.5*t);
end

function [bcx, bcy] = boundary(bcx,bcy,x,t)
m = length(x);

bcx(1:m)=exp(-x(1)-x-0.5*t); bcx(m+1:2*m)=exp(x(1)+x+0.5*t);
bcx(2*m+1:3*m)=exp(-x(end)-x-0.5*t); bcx(3*m+1:end)=exp(x(end)+x+0.5*t);

bcy(1:m)=exp(-x-x(1)-0.5*t); bcy(m+1:2*m)=exp(x+x(1)+0.5*t);
bcy(2*m+1:3*m)=exp(-x-x(end)-0.5*t); bcy(3*m+1:end)=exp(x+x(end)+0.5*t);

end

function u2v = U2V(W,n)
u2v = (W(1:n).^2).*W(n+1:2*n);
u2v = [u2v; u2v];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Better to work with t = 20


% Video stuff for this file
%{
%________________________________________
video_on = 0;
%________________________________________
if video_on && (t>=0)&&(mod(freq,100)==0)
    
       % surf(X,Y,reshape(V(1:n),m,m));  % <<---- u
        surf(X,Y,reshape(V(n+1:2*n),m,m));   % <<---- v
      %  surf(X,Y,reshape(u_exact(X,Y,t),m,m) )
        shading interp
       % view(90,0)
        colorbar
        title(['Numerical solution at t = ',num2str(t)]);
        %theAxes=[x_l x_r y_d y_u];
        axis(theAxes);
        grid;xlabel('x');
        legend('u')
        ax = gca;          % current axes
        ax.FontSize = 16;
        currFrame = getframe;
       % writeVideo(vidObj,currFrame);
    end

%________________________________________
if video_on
    theAxes=[x_l x_r x_l x_r 0 15]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
   % vidObj = VideoWriter('System');
    open(vidObj);
end
%________________________________________
if video_on
    close(vidObj);
end
%________________________________________
%}
