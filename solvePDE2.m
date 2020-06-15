function [t, C] = solvePDE(f, v, D, tRange, xc, ~)

D=1;
solvePDE(@(t,x,C) , @(t,x) 0; @(t,x) 1; [0 10]; linspace(0,10,100), exp(-((1:100)-50).^2);
%Solution of an advection-reaction-diffusion equation:
xc = linspace(0,10,100);
[t, C] = solvePDE(@(t,x,C) 0.01*C, @(t,x) 0, @(t,x) 1, [0 100], xc , exp(-(xc-5).^2));
clf
surface(t,xc,C')
shading interp


if nargout==0
    clf
    surf(t,xc,C')
    xlabel('t')
    ylabel('x')
    zlabel('C')
    shading interp
end

    function dCdt = derivative(t,C)
        
        % Initialize
        dCdt = zeros(n,1);  % Rate of change of average concentration in each cells
        
        %
        % Transport: Loop over all interfaces
        %
        for i = 1:(n-1)
            %
            % Compute flux from cell i to cell i+1
            %
            %
            % Advective flux: Upwind method
            vi = v(t,xi(i));
            
            if vi>0
                Ja = vi*C(i);
            else
                Ja = vi*C(i+1);
            end
            
            Jd = D(t,xi(i))*(C(i)-C(i+1))/dx;   % diffusive flux
            J = Ja  + Jd;         % Total flux
            %
            % Remove transport from cell i and add to cell i+1
            %
            dCdt(i) = dCdt(i) - J/dx;
            dCdt(i+1) = dCdt(i+1) + J/dx;
        end        
        %
        % Reaction: Loop over all cells
        %
        for i = 1:n
            dCdt(i) = dCdt(i) + f(t,xc(i),C(i));
        end
    end
end
