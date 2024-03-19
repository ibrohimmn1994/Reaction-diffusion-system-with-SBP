%{
Finite element analysis of nonlinear reaction–diffusion system of
Fitzhugh–Nagumo type with Robin boundary conditions



%}
m = 201;
x_l=-20.0;x_r=20.0;
len = x_r - x_l;                  
n=m*m;
n_t = 2*n;
h=(x_r-x_l)/(m-1);


CFL=0.0005;
k=CFL*h^2

K = 2;
Ik=speye(K);
Im=speye(m);
In=speye(n);
I=speye(K*m*m);
I2 = eye(2);

e1=[1 0];
e2=[0 1];

SBP4;

Dxx = kron(sparse(D2), Im); Dyy = kron(Im, sparse(D2));
HIx = kron(sparse(HI), Im); HIy = kron(Im, sparse(HI));
HI = HIx*HIy;
HI_bar = sparse(kron(Ik, HI)) ;

Hx = kron(sparse(H), Im); Hy = kron(Im, sparse(H));
H_bar = Hx*Hy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear HIx HIy D1 D2 M Q M_U Q_U S_U H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_E = kron(sparse(d_m), Im); d_W = kron(sparse(d_1), Im);
d_N = kron(Im, sparse(d_m)); d_S = kron(Im, sparse(d_1));


Lx = [ kron(I2,d_W); kron(I2,d_E)];
Ly = [kron(I2,d_S) ; kron(I2,d_N)];

Px = I - HI_bar*Lx'*((Lx*HI_bar*Lx')\Lx);
Py = I - HI_bar*Ly'*((Ly*HI_bar*Ly')\Ly);
P = Px*Py;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxx = kron(Ik,Dxx);
Dyy = kron(Ik,Dyy);


C1 = ones(n_t,1); C1(n+1:n_t) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=linspace(x_l,x_r,m);	% discrete x values
[X,Y] = meshgrid(x,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

video_on = 1;

t_1= 5;
max_itter=floor(t_1/k);
V = zeros(2*n,1); 

V(1:n) =  reshape(u_initial(X,Y),n,1);
V(n+1:2*n) = reshape(v_initial(X,Y),n,1);
%__________________________________________________________________________
if video_on
    
    theAxes=[x_l x_r x_l x_r -1 1.5]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
   % vidObj = VideoWriter('System');
    open(vidObj);
end


t=0.0;
freq = 0;

for nr_itter=1:max_itter
    % stage1
    W = V ;
  
    w1 =   C1.*(Dxx*W + Dyy*W) +  F(W,n);
    %______________________________________________________________________
    %stage2
    

    W = P*(V+k/2*w1) ;
    w2 = C1.*(Dxx*W + Dyy*W) + F(W,n) ; 
    %______________________________________________________________________
    %stage3
    
   
    W = P*(V+k/2*w2) ;
    w3 =  C1.*(Dxx*W + Dyy*W) + F(W,n);  
   
    %______________________________________________________________________
    %stage 4
    

    W = P*(V+k*w3) ;
    w4 =  C1.*(Dxx*W + Dyy*W) + F(W,n) ;  
    %_____________________________________________________________________
    V= P*(V+k/6*(w1+2*w2+2*w3+w4)) ;
    
    
    t=t+k;

    

    if video_on && (t>=0)&&(mod(freq,500)==0)
        t
        surf(X,Y,reshape(V(1:n),m,m));  % <<---- u
       % surf(X,Y,reshape(V(n+1:2*n),m,m));   % <<---- v
      %  surf(X,Y,reshape(u_exact(X,Y,t),m,m) )
       
        shading interp
        view(0,90)
        colorbar
       % view(90,0)
      
        title(['Numerical solution at t = ',num2str(t)]);
        
        axis(theAxes);
        grid;xlabel('x');
        legend('u')
        ax = gca;          % current axes
        ax.FontSize = 16;
        currFrame = getframe;
       % writeVideo(vidObj,currFrame);
    end
    freq = freq + 1;
end


function vi = v_initial(x,y)
for i=1:length(x)
    for j=1:length(y)
        if (x(i,j)<1 && y(i,j)<-10)
            vi(i,j) = 0.15;
        else
            vi(i,j) = 0.0;
        end
    end
end    
end

function ui = u_initial(x,y)
for i=1:length(x)
    for j=1:length(y)
        if (x(i,j)<0 || y(i,j)>5)
            ui(i,j) = 0.0;
        else
            ui(i,j) = (1+exp(4*(abs(x(i,j))-5)))^(-2) - (1+exp(4*(abs(x(i,j))-1)))^(-2) ;
        end
    end
end  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = F(W,n)
gamma = 0.5;alpha=1; a=0.1; eps =0.001; 
u = W(1:n);
v = W(n+1:2*n);
f = 1/eps*(-u.*(u-1).*(u-a) - v);
g = alpha*u - gamma*v;
R = [f;g];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Better to work with t = 20
