

function  DE = discrete_equation(U, DE_coeff, name,x,t)

switch name
    
    case "solitons"
        [su, sv] = MMS.source(x,t);
        DE = DE_coeff*(U) + [SF_2.F(U);SF_2.G(U)] + [su;sv];

    case "analytical_2"
        m = length(U)/2;
        vu2 =  U(1:m).*(U(m+1:2*m)).^2;   
        DE = DE_coeff{1}*(U) + DE_coeff{2}*([vu2;vu2]) + DE_coeff{3}*(U) + DE_coeff{4};
    case "FHN"
        n = length(U)/2;
        u = U(1:n);
        v = U(n+1:2*n);
        f = -u.*(u-1).*(u-0.15) - v;
        g = 0.005*u - 0.025*v;
        R = [f;g];
        S = MMS.sourse_FHN(x{1},x{2},t);
        DE = DE_coeff{1}.*(DE_coeff{2}*U + DE_coeff{3}*U) +  R + S;

    case "Brusselator"
        n = length(U)/2;
    
        u = U(1:n);
        v = U(n+1:2*n);
        u2v = (u.^2).*v;
        u2v = [u2v; u2v];
        
        DE = DE_coeff{1}*(DE_coeff{2}*U + DE_coeff{3}*U) + DE_coeff{4}.*u2v + DE_coeff{5}.*[u;u];

    
    otherwise
        return
end

end





    