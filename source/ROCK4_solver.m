function [V_new,t_now ]= ROCK4_solver(t_begin, t_end,coeff,BC_coeff, V, P , BC, k, name,x,atol,rtol)


spectralRadius = @SpectralRadius;

%Rock_step = @rock_step;
Rock_step = @rock_step_analytical2;
%Rock_step = @rock_step_with_boundary;
%atol = 1e-6;%1e-4;
%rtol = 1e-6;%1e-3;

t_step = k;
t_now = t_begin;

safety = 0.9;%0.8
%scale = 1.0;

err_loc=0.0;
err_loc_old=0.0;

%temp1=1.0;
%temp2=1.0;

nrej=0;
naccp=0;
nrho=0;
%fac = 1;
% facmax = 5;
facmax = 1.5;
%mdeg = 5;
rho = 0.0;

%spectrOK = False;
EPS = 1.0e-16;
hmin = EPS * 1.0e1 * max(t_begin, t_end);

t_step = min(t_step, t_end-t_now);

rho = spectralRadius(  V, rho,coeff, name,x,t_begin); %<<-----------------------

t_step = max(0.25./rho, t_step);
t_step = min(t_step, t_end-t_now);
t_step_old = t_step;


while (t_now < t_end)
    %disp(t_step)
   % disp(mdeg)
    
    if (mod(nrho, 5)== 0)
        rho = spectralRadius(  V, rho,coeff, name,x,t_now);%<<----------------------
    end
    nrho = nrho +1;
    mdeg = 1 + round( sqrt( (3+t_step*rho)./0.353) );
   
    if (mdeg>152)
        t_step = safety * (152^2 * 0.353 - 3)./rho;
        mdeg = 152;
    end
    mdeg = max(mdeg, 5)-4;
    t_step = min(t_step, t_end-t_now);
    %disp(t_step)
  
    
    [V_new, err_loc] = Rock_step( mdeg, atol, rtol, coeff, V, P , BC, t_step, name,x,t_now); %<<-----------------
   % [V_new, err_loc]= Rock_step( mdeg, atol, rtol, coeff,BC_coeff, V, P ,BC, t_step, name, x, t_now);% this is for the analytical

    
    
   % disp(t_step)
    %disp(err_loc)
    %{
    m=101;
    x=linspace(0,1,m)';	
    figure(i);
    plot(x, V_new(1:m))
    %}

    if (err_loc >= 1)%0.43
        %disp('rej')
        nrej  = nrej + 1;
        t_step = safety *t_step./(err_loc.^(1/5));
        t_step = max(t_step, hmin);

    else
   
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %disp(t_step)
        naccp = naccp + 1;
        t_now = t_now + t_step;%<<--------------------
        V = V_new;

        fac = facmax;
        if (naccp == 1)
            
            temp2 = err_loc.^(2.0/5.0);
            if (safety <fac*temp2)
                fac = safety./temp2;
            end
       
        else
            temp1 = safety *t_step * err_loc_old.^(1.0/5.0);
            temp2 = t_step_old * err_loc.^(2.0/5.0);
            if (temp1 < fac* temp2)
                fac = temp1./temp2;
            end

        end
        
        t_step_old = t_step;
        err_loc_old = err_loc;
        scale = max(0.1, fac);
        t_step = t_step *scale;
        %disp(scale)
        t_step = max(hmin, t_step);
    end
end
%disp(nrej)
%disp(naccp)



end