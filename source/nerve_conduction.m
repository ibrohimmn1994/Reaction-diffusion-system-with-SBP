video_on = 1;


global m
m = 399;
x_l=0;x_r=100;
len = x_r - x_l;                    % Number of gridpoints
n=2*m;
h=(x_r-x_l)/(m-1);
%__________________________________________________________________________
if video_on
    theAxes=[x_l x_r -1 1]; 
    %theAxes=[x_l x_r 0.0019 8];%.0021]; 
    %theAxes=[x_l x_r 0 1]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
    vidObj = VideoWriter('System');
    open(vidObj);
end




t_1= 200;                    % End time
CFL=0.1;
k=CFL*h^2;

I2=speye(2);
I=speye(n);
e1=[1 0];
e2=[0 1];

SBP4;


HII = kron(I2,HI);
L = [ kron(e1,d_1); kron(e1,d_m)];

gg = zeros(2,1);
gg(1) = -0.3;
P = I - HII*L'*((L*HII*L')\L);
BC = HII*L'*((L*HII*L')\gg); 

eta = 0.008;
beta = 2.54;

% the terms
T1 = sparse([D2 0*D2;...
    0*D2 0*D2]);
T2 = sparse(diag([-ones(m,1); ...
    eta*ones(m,1)]) );
T3 = sparse(diag([-ones(m,1);...
    -eta*beta*ones(m,1)]));

coeff = {T1,T2,T3};
max_itter=floor(t_1/k);
x=linspace(x_l,x_r,m)';	
V=sparse(zeros(n,1) );          



t=0.0;
freq = 0;
for nr_itter=1:max_itter
    % For ROCK4 use this line
    %MyRock4_nerve_conduction;
      
   % V = MyRock4(coeff, V, P , BC, 152, k,"nerve_conduction");

    % For RK4 use this block
    
    U = P*V + BC;
    w1 = P*(T1*(U) + T2*[f(U(1:m));U(1:m)] + T3*[U(m+1:2*m);U(m+1:2*m)] );
    U = P*(V+k/2*w1)+ BC;
    w2 = P*(T1*(U) + T2*[f(U(1:m));U(1:m)] + T3*[U(m+1:2*m);U(m+1:2*m)] );
    U = P*(V+k/2*w2)+ BC;
    w3 =  P*(T1*(U) + T2*[f(U(1:m));U(1:m)] + T3*[U(m+1:2*m);U(m+1:2*m)] );
    U = P*(V+k*w3)  + BC;
    w4 =  P*(T1*(U) + T2*[f(U(1:m));U(1:m)] + T3*[U(m+1:2*m);U(m+1:2*m)] );
    V=V+k/6*(w1+2*w2+2*w3+w4);
    
    
    t=t+k;

    if video_on && (t>=0.0)&&(mod(freq,1000)==0)
       % freq

        plot(x,real(V(1:m)),'r',x,real(V(m+1:n)),'b--','LineWidth',1);
        %title(['Numerical solution at t = ',num2str(t)]);
        %theAxes=[x_l x_r y_d y_u];
        axis(theAxes);
        grid;xlabel('x');
        %legend('v^{(1)}','v^{(2)}')
        ax = gca;          % current axes
        ax.FontSize = 16;
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    freq = freq + 1;
end

if video_on
    close(vidObj);
end
%}
function F =  f(u)
alpha = 0.139;
F = u.*(u-alpha).*(u-1);
end
