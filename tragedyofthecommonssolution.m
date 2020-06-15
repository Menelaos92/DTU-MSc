%
% Illustration of the tragedy of the commons
%
% The resource has logistic growth:
%  N' = r N(1-N/K) - E N
%
% The profit is the revenue minus the cost:
%  R = p*E*N-c*E
%
% The increase in effort is proportional to the profit per unit effort (R/E):
%
%  E' = k*R/E = k*(p*N-c)
%
% At equilibrium this gives Eequ = r*(1-c/(p*k)) and Nequ = c/p.
%
% What is the revenue?
% What is the dynamics when k is slower (faster) than r?
%


%
% Define parameter values:
%
K = 1000; % Carrying capacity
r = 0.1;  % Population growth rate
p = 15;   % Price ($ per capita)
c = 10000; %.2 * p*K/2; % Cost per effort. 
                     % If we assume that the revenue is large at Emsy, then c << p*K/2 
k = 1 * r^2/(K*p-c); % Speed of change of effort. 
                     % If the constant in front is =1, then the effort changes
                     % on the same time scale as the population dynamics.
%
% Solve the ODEs:
% Start at N=K (carrying capacity) and no effort.
%
clf
[t y] = ode45(@derivative, [0 150/r], [K, 0],[],K,r,p,c,k);
%
% Define temporary variables to increase readability
%
N = y(:,1);
E = y(:,2);
%
% Plot the results
%
onez = ones(length(t),1); % Make a vector of ones
Eequ = r*(1-c/(p*K));     % Equilibrium solution of effort
Nequ = c/p;               % Equilibrium solution of abundance

clf
% Plot Abundance
subplot(3,1,1)
plot(t,N,'-', t,Nequ.*onez,'k-', t,K/2.*onez,'r--', 'linewidth',2)
title('Abundance')
legend('Solution','Equ. solution','MSY solution')

% Plot profit per effort:
subplot(3,1,2)
plot(t,p*N-c,'b-', t,0*t,'k-', t,onez*p*K/2-c*r/2,'r--', 'linewidth',2)
title('Revenue per unit effort (R/E)')

% Plot effort:
subplot(3,1,3)
plot(t,E,'-',t,onez*Eequ,'k-', t,r/2.*onez,'r--', 'linewidth',2)
title('Effort')
xlabel('Time')


    function dydt = derivative(t,y,K,r,p,c,k)
        % Define temporary variables to increase readability
        N = y(1);
        E = y(2);
        % Calculate the derivatives:
        dNdt = r*(1-N/K)*N - E*N;
        dEdt = k*(p*N-c);
        % Make sure that the effort cannot become negative:
        if E<=0 && dEdt<0
          dEdt = 0; 
        end
        % Assemble the vector with the derivatives:
        dydt = [dNdt dEdt]';
    end