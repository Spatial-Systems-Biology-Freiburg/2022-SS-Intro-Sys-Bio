function TS = lsa(k)
% linear stability analysis
dpq = FourierModes(20,20);
TS = 0;
% evaluate single cell (non-spatial) stability
% solve ODE:
[t,y] = ode15s(@odes, [0 10000], zeros(5,1)+0.0001, [], 1, 0, k);
ss = y(end,:); % get steady state solutions
% Assign steady state values:
Ac = ss(1); An = ss(2); Ic = ss(3); In = ss(4); T = ss(5);

% Fill in Jacobian with steady state solution and parameters (k):
J = [-k(1)-k(3),k(1) + T*k(2),0,0,An*k(2);
  k(1), (2*An*k(5))/(In + k(6))-k(3)-T*k(2)-k(1),0,-(An^2*k(5))/(In + k(6))^2, -An*k(2);
  0,0,-k(1)-k(7)-T*k(8),k(1),-Ic*k(8);
  0,2*An*k(10),k(1) + T*k(8),-k(1)-k(7),Ic*k(8);
  0,-T*k(2),-T*k(8),0, -k(12)-An*k(2)-Ic*k(8)];

D = diag([1 0 k(9) 0 0]); % Diagonal matrix with diffusion coefficients

%Find positive eigenvalues for J:
pos_eig = find(eig(J)>0);
% Non-spatial system should be stable, so reject if positive eigenvalues:
if (numel(pos_eig) >= 1)
    return
end

%% Check stability for spatial problem
for idx = length(dpq):-1:1
    % Jacobian - diffusion matrix * wavemodes:
    Apq = J-D*dpq(idx);
    % Find positive eigenvalues:
    evals = eig(Apq);
    pos_eig = find(evals>0);
    % If for any given p,q the matrix Apq has an eigenvalue
    % with a positive real part the system is unstable:
    if (numel(pos_eig) >= 1)
        TS = 1;
        break;
    end
end