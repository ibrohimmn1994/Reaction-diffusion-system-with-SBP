function rho = SpectralRadius( yn, rho ,coeff,  name, x, t)

F = @discrete_equation;
len = length(yn);
%converged= 0;
lambda = 0;
itmax = 30;%should be 30
nfune = 0;
EPS = 3.0e-16;


%funn = P*F(yn,coeff, name, x, t); %fz  <<---------------------
funn = F(yn,coeff, name, x, t); %fz  <<---------------------

nfune = nfune + 1;

v= funn; %z

vnrm = sqrt(sum(v.*v));%znor
ynrm = sqrt(sum(yn.*yn));
%disp('-----')
if (abs(ynrm)>0 && abs(vnrm) >0)
   % disp('1')
    dynrm = ynrm *sqrt(EPS);
    v = yn + (dynrm/vnrm)*v;

else
    %disp(abs(ynrm));
    %disp(abs(vnrm));

    if (abs(ynrm)>0)
       % disp('2')
        dynrm = ynrm *sqrt(EPS);
        v = (1 + sqrt(EPS))*yn;
    end
    if (abs(vnrm) > 0)
       % disp('3')
        dynrm = EPS;
        v = (dynrm/vnrm) * v;
    end
    if (abs(ynrm)==0 && abs(vnrm) == 0)
       % disp('4')
        dynrm = EPS;
        v(1) = dynrm; % <<<<----
    end
end


for iter=1:itmax
   % funv = P*F(v,coeff, name,x,t); %<<-------------------------
    funv = F(v,coeff, name,x,t); %<<-------------------------
    nfune = nfune + 1;
    dfunnrm = sqrt(sum((funv - funn).*(funv - funn)) );
    lambda_old = lambda;
    lambda = dfunnrm./dynrm;
    %safety instead of 1.1
    rho = 1.1*lambda;
    % if iter >=2 and <= 0.05*rho or lambda
    if ( abs(lambda-lambda_old)<=(0.01*max(lambda,lambda_old)))
       % disp('a')
        % this we do not need it, we can compute the perturbation as above
        % or alternatievly we can store ev and use it as a perturbation
        % in the next call of the implicit spectral method
 %       ev = v - yn; % <<<----------------------------- check again here
     %   converged = 1;
   %     break;
        return;
    end
    if ( abs(dfunnrm) > 0 )
        v = yn + (funv - funn) * (dynrm/dfunnrm);
    else
        index = mod(iter,len );
        v(index) = yn(index) - (v(index) - yn(index));
    end



end

end