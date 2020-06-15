% =========================================================================
% solvePDEs Solve a system of advection-diffusion-reaction equations numerically.
% 
% [tv xc C] = solvePDEs(v,D,f,xi,tv,C0,Fickian)
% 
% solves the PDEs specified by v, D, f, for t in tv, with initial condition C0.
%        
% The numerical method is semi-implicit finite volume/finite difference.
%
% Input:    v : A function v(t,x,C,dCdx) giving advection
%           D : A function D(t,x,C) giving diffusion
%           f : A function f(t,x,C) giving reaction
%           xi: A vector of interfaces between grid cells
%           tv: A vector of time points
%           C0: A function C0(x) giving initial condition
%           Fickian: 1 for Fickian diffusion, 0 for Random Walk diffusion
%
%           All functions should return vector outputs, one for each PDE.
%
% Output:   tv: vector of time points
%           xc: vector of spatial grid points
%           C: the solution, a 3D array, with time, space, and component.
%
% Examples:
% 
% solvePDEs(@(x,t,C,dCdx)1,@(t,x,C)1,@(t,x,C)0,0:0.05:1,0:0.01:0.5,@(x)1+cos(pi*x));
%
% solvePDEs(@(x,t,C,dCdx)[0,0],@(t,x,C)[1,1],@(t,x,C)C.*(1-C),0:50,0:0.5:50,@(x)[double(x<1),double(x>49)]);

% ========================================================================
function [tv xc C] = solvePDEs(v,D,f,xi,tv,C0,Fickian)

if nargin<7,
    Fickian = 1;
end

% ----------------------------------------------
%  Numerical grid:
% ----------------------------------------------
n = length(xi)-1;         % no. of grid cells
if n == 1,
    n = 51;
    xi = linspace(xi(1),xi(2),101);
end

m = size(C0(xi(1)),2);

nt = length(tv);
dt = tv(2)-tv(1);

dx = diff(xi);
dx = dx(1);
xc = xi(2:end)-0.5*dx; % Center points of grid cells on the x-axis

% ----------------------------------------------
% Solve the system of ODEs:
% ----------------------------------------------
tic

C = zeros(nt,n,m);

for i=1:n,
    C(1,i,:) = C0(xc(i));
end

for i=2:nt,
    % Implicit advection-diffusion step
    for j=1:m,
        % Setup I - G'*dt
        IGt = sparse(1:n,1:n,ones(1,n),n,n,3*n);
        
        % Loop over all interfaces
        for k=1:(n-1),
            if Fickian,
                % Evaluate diffusivity at the interface
                DD = D(tv(i),xi(k+1),0.5*(C(i-1,k,:)+C(i-1,k+1,:))) / dx / dx * dt;
                DD = DD(j);
                
                IGt(k,k) = IGt(k,k) + DD;
                IGt(k+1,k) = IGt(k+1,k) - DD;
                IGt(k,k+1) = IGt(k,k+1) - DD;
                IGt(k+1,k+1) = IGt(k+1,k+1) + DD;
            else
                % Evaluate diffusivity at the left cell             
                DD = D(tv(i),xc(k),C(i-1,k,:)) / dx / dx * dt;
                DD = DD(j);
                
                IGt(k,k) = IGt(k,k) + DD;
                IGt(k+1,k) = IGt(k+1,k) - DD;
                
                % Evaluate diffusivity at the right cell             
                DD = D(tv(i),xc(k+1),C(i-1,k+1,:)) / dx / dx * dt;
                DD = DD(j);
                
                IGt(k,k+1) = IGt(k,k+1) - DD;
                IGt(k+1,k+1) = IGt(k+1,k+1) + DD;
            end
            
            % Evaluate velocity at the interface
            vv = v(tv(i-1),xi(k+1),0.5*(C(i-1,k,:)+C(i-1,k+1,:)),(C(i-1,k,:)-C(i-1,k+1,:))/dx) / dx * dt;
            vv = vv(j);
            
            % Upwind method
            if vv>0,
                IGt(k,k) = IGt(k,k) + vv;
                IGt(k+1,k) = IGt(k+1,k) - vv;
            else
                IGt(k+1,k+1) = IGt(k+1,k+1) - vv;
                IGt(k,k+1) = IGt(k,k+1) + vv;
            end
        end
        
        % Implicit step solving transport
        C(i,:,j) = IGt \ C(i-1,:,j)';
    end
    
    % Predictor-corrector step for reaction dynamics
    for k=1:n,
        % Reaction dynamics evaluated at the old step 
        f1 = squeeze(f(tv(i-1),xc(k),squeeze(C(i-1,k,:))'));
        f1 = f1(:);
        
        % Predict new density, using both transport and reaction
        Cp0 = squeeze(C(i,k,:));
        Cp0 = Cp0(:);
        Cp = Cp0 + f1*dt;
        
        % Evaluate dynamics evaluated at the prediction
        f2 = squeeze(f(tv(i),xc(k),Cp'));
        f2 = f2(:);
        
        % Use a compromize between the two reaction dynamics
        Cp = Cp0 + 0.5*(f1+f2)*dt;
        C(i,k,:) = Cp;
    end
    
end
    

% subplot(2,2,1)
% plot(C0(1:n))
% subplot(2,2,2)
% plot(C0(n+(1:n)))
% dCdt = derivative(0,C0(:));
% subplot(2,2,3)
% plot(dCdt(1:n))
% subplot(2,2,4)
% plot(dCdt(n+(1:n)))

% [t C] = ode45(@derivative, tspan, C0(:));
% toc

% C = reshape(C,[length(t),n,m]);

% Extract solution fields

if nargout==0,
    for i=1:m
        subplot(1,m,i)
        imagesc(tv,xc,C(:,:,i)')
        xlabel('t')
        ylabel('x')
        title(i)
    end
end
