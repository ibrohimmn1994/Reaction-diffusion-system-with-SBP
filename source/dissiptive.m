video_on = 1;

t_1= 100;
max_itter=floor(t_1/k);

UU = reshape(initial_UU(X,Y),n,1);


%rings
%{
tau = 4.1; q=0.1; Du = 5.3e-; Dv = 0;Dw = 9.64e-3;
lambda = 0.95; k1 = -0.08; k3 = 0.25; k4 =1.0;
%}

% rings with ID

% 2 echelons

tau = 48.1; q=1.0; Du = 5.3e-4; Dv = 0;Dw = 9.64e-3;
lambda = 0.95; k1 = -0.08; k3 = 0.25; k4 =1.0;

% one echelons
%{
tau = 48; q=0.48; Du = 5.3e-5; Dv = 0.0;Dw = 9.64e-2;
lambda = 0.95; k1 = -0.08; k3 = 0.25; k4 =1.0;
%}


% breathing
%{
tau = 1.1; q=0.5; Du = 1e-3; Dv = 0;Dw = 0.01;
lambda = 5.67; k1 = -1.04; k3 = 1; k4 = 3.33;
%}


VV = ones(n,1)*v0;
WW = ones(n,1)*w0;


%__________________________________________________________________________
if video_on
    theAxes=[x_l x_r x_l x_r -3 6];
    %theAxes=[x_l x_r 0.0 30]; 
    %theAxes=[x_l x_r 0.0019 8];%.0021]; 
    %theAxes=[x_l x_r 0 1]; 
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 0 scrsz(3)/2 scrsz(4)-1])
    vidObj = VideoWriter('System');
    open(vidObj);
end
%P=1;gg=0;


t=0.0;
freq = 0;
%DY=0;

BCx_u = 0;BCx_v = 0;BCx_w = 0;
Px_u =1; Px_v =1; Px_w =1;
Py_u =1; Py_v =1; Py_w =1;

for nr_itter=1:max_itter
    UUx = Px_u*UU+BCx_u; VVx = Px_v*VV+BCx_v; WWx =Px_w*WW+BCx_w; 
    UUy = Py_u*UU+BCy_u; VVy = Py_v*VV+BCy_v; WWy =Py_w*WW+BCy_w; 
  
    w1_u =  Px_u*( Du*DX*(UUx) )  + Py_u*( Du*DY*(UUy)) +lambda*UU -UU.^3 -k3*VV -k4*WW + k1;
    w1_v =  Px_v*( 1/tau*Dv*DX*(VVx) )  + Py_v*( 1/tau*Dv*DY*(VVy)) +1/tau*UU -1/tau*VV;
    w1_w =  Px_w*( 1/q*Dw*DX*(WWx) )  + Py_w*( 1/q*Dw*DY*(WWy)) +1/q*UU - 1/q*WW;

    UU0 = UU+k/2*w1_u; VV0 = VV+k/2*w1_v; WW0 = WW+k/2*w1_w; 
    UUx = Px_u*(UU0)+BCx_u; VVx = Px_v*(VV0)+BCx_v; WWx =Px_w*(WW0)+BCx_w; 
    UUy = Py_u*(UU0)+BCy_u; VVy = Py_v*(VV0)+BCy_v; WWy =Py_w*(WW0)+BCy_w; 
    

    w2_u =  Px_u*( Du*DX*(UUx) )  + Py_u*( Du*DY*(UUy)) +lambda*UU0 - UU0.^3 -k3*VV0 -k4*WW0 + k1;
    w2_v =  Px_v*( 1/tau*Dv*DX*(VVx) )  + Py_v*( 1/tau*Dv*DY*(VVy)) -1/tau*UU0 -1/tau*VV0;
    w2_w =  Px_w*( 1/q*Dw*DX*(WWx) )  + Py_w*( 1/q*Dw*DY*(WWy)) +1/q*UU0 - 1/q*WW0;

    UU0 = UU+k/2*w2_u; VV0 = VV+k/2*w2_v; WW0 = WW+k/2*w2_w;
    UUx = Px_u*(UU0)+BCx_u; VVx = Px_v*(VV0)+BCx_v; WWx =Px_w*(WW0)+BCx_w; 
    UUy = Py_u*(UU0)+BCy_u; VVy = Py_v*(VV0)+BCy_v; WWy =Py_w*(WW0)+BCy_w; 

    w3_u =  Px_u*( Du*DX*(UUx) )  + Py_u*( Du*DY*(UUy))+lambda*UU0 -UU0.^3 -k3*VV0 -k4*WW0 + k1;
    w3_v =  Px_v*( 1/tau*Dv*DX*(VVx) )  + Py_v*( 1/tau*Dv*DY*(VVy)) +1/tau*UU0 -1/tau*VV0;
    w3_w =  Px_w*( 1/q*Dw*DX*(WWx) )  + Py_w*( 1/q*Dw*DY*(WWy)) +1/q*UU0 - 1/q*WW0;


    UU0 = UU+k*w3_u; VV0 = VV+k*w3_v; WW0 = WW+k*w3_w;
    UUx = Px_u*(UU0)+BCx_u; VVx = Px_v*(VV0)+BCx_v; WWx =Px_w*(WW0)+BCx_w; 
    UUy = Py_u*(UU0)+BCy_u; VVy = Py_v*(VV0)+BCy_v; WWy =Py_w*(WW0)+BCy_w; 
    
    w4_u =  Px_u*( Du*DX*(UUx) )  + Py_u*( Du*DY*(UUy))+lambda*UU0 -UU0.^3 -k3*VV0 -k4*WW0 + k1;
    w4_v =  Px_v*( 1/tau*Dv*DX*(VVx) )  + Py_v*( 1/tau*Dv*DY*(VVy)) +1/tau*UU0 -1/tau*VV0;
    w4_w =  Px_w*( 1/q*Dw*DX*(WWx) )  + Py_w*( 1/q*Dw*DY*(WWy)) +1/q*UU0 - 1/q*WW0; 

  
    UU=UU+k/6*(w1_u+2*w2_u+2*w3_u+w4_u);
    VV=VV+k/6*(w1_v+2*w2_v+2*w3_v+w4_v) ;
    WW=WW+k/6*(w1_w+2*w2_w+2*w3_w+w4_w) ;
    
   
    
    t=t+k;

    if video_on && (t>=0.0)&&(mod(freq,2000)==0)
        %freq
       % V(n+1:2*n)
        surf(X,Y,real(reshape(UU,m,m)));%+reshape(V(n+1:2*n),m,m)
       % hold on
       
        
        shading interp
        view(0,90)
        colorbar
        title([' t = ',num2str(t)]);
        %theAxes=[x_l x_r y_d y_u];
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

if video_on
    close(vidObj);
end

%__________________________________________________________________________
%{
figure(2);
plot(x,V(1:m), x,V(m+1:n),'LineWidth',1);
%plot(x,V(1:m,1),'r','LineWidth',1);
title(['Numerical solution at t = ',num2str(t_1)]);
legend('V^{1}','V^{2}')
xlabel('x');ylabel('V');
ax = gca;          % current axes
ax.FontSize = 16;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function ID = initial_VV(x,y)
global v0
%u0 = -1.5;
ID = v0 + 3*exp( -((x-0.02)/0.05).^2 -((y)/0.05).^2);
end
function ID = initial_WW(x,y)
global w0
%u0 = -1.5;
ID = w0 +  ( 2.5*exp( -((x-0.1)/0.05).^2-((y)/0.05).^2) ) ;
end

function ID = initial_UU(x,y)
global u0
ID = u0 +   ( 200*exp( -((x)/0.05).^2-((y)/0.05).^2) ) - ( 200*exp( -((x+0.01)/0.05).^2-((y)/0.05).^2) ) ;
end

%}