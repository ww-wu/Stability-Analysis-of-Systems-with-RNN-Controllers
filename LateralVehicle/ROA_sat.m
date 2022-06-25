function [Px, Px_trace] = ROA_sat(sys,Whx,Whh,P_hat)

A_G = sys.A;
B_G = sys.B;
n_G = size(A_G,1);
C_G = sys.C;
n_u = size(B_G,2);
n_y = size(C_G,1);

L = length(Whh);
n_h = zeros(L,1);
for e = 1:L
    n_h(e) = size(Whx{e},1);
end
n = sum(n_h);

Nuh = [zeros(n_u, n-n_h(L)) , Whx{L+1}];
Nvy = [Whx{1} ; zeros(n-n_h(1), n_y)];
Nvh = [zeros(n_h(1),n) ; blkdiag(Whx{2:L}) , zeros(n-n_h(1), n_h(L))];
Nvd = blkdiag(Whh{:});

[v_max,v_min,h_max,h_min,alpha_max,alpha_min,sigma_max,sigma_min] = interval_bound_propagation(C_G, Whx, Whh, P_hat);

% IQC for saturation
u_max = abs(Whx{L+1}) * h_max{L};
umax = 30/180*pi;
sec_u = min(1 , umax/u_max);
slop_u = u_max < umax;
Bsat1 = -1;
Bsat2 = 1;
Csat = [0 ; 0 ; 1 ; 0];
Dsat1 = [1 ; -sec_u ; 1 ; -slop_u];
Dsat2 = [-1 ; 1 ; -1 ; 1];
% The state is sat(u)-u
sat_bound_square = max(u_max-umax, 0)^2;

% IQC for LTI uncertainty
delta_norm = 0.1;
dt = 0.02;
A_Delta = -eye(2);
n_Delta = size(A_Delta,1);
A_Delta = A_Delta*dt + eye(n_Delta);
B_Delta = eye(2);
B_Delta = B_Delta*dt;
C_Delta = [0,1,0,0;0,0,0,1]';
D_Delta = [[1,0,0,0]' , [0,0,1,0]'];
% The upper bound of the state is the upper bound of its inputs,
% which are the outputs of the saturation and LTI uncertainty respectively.
delta_bound_square = (min(u_max, umax) * [1 ; delta_norm]).^2;

v_max = vertcat(v_max{:});
v_min = vertcat(v_min{:});
h_max = vertcat(h_max{:});
h_min = vertcat(h_min{:});
v_bound_square = max(v_max.^2, v_min.^2);
h_bound_square = max(h_max.^2, h_min.^2);
alpha_max = horzcat(alpha_max{:})';
alpha_min = horzcat(alpha_min{:})';
sigma_max = horzcat(sigma_max{:})';
sigma_min = horzcat(sigma_min{:})';

B_Phi_v = -diag(sigma_max);
B_Phi_h = eye(n);
C_Phi = [zeros(2*n, n) ; eye(n) ; zeros(n)];
D_Phi_v = [diag(alpha_max) ; -diag(alpha_min);...
           diag(sigma_max) ; -diag(sigma_min)];
D_Phi_h = [-eye(n) ; eye(n) ; -eye(n) ; eye(n)];

% combination
n_Xi = n_G + 2*n + n_Delta + 1;
A_Xi = blkdiag([A_G , zeros(n_G,2*n) ; Nvy*C_G , zeros(n,n) , Nvd ; zeros(n,n_G+n+n)], A_Delta, 0);
B_Xi = [zeros(n_G,n) , B_G , B_G ;...
        Nvh , zeros(n,2) ;...
        eye(n) ,zeros(n,2) ;...
        zeros(n_Delta, n) , B_Delta ;...
        Bsat1*Nuh , Bsat2 , 0];
C_Xi = blkdiag([D_Phi_v*Nvy*C_G , C_Phi*B_Phi_v , D_Phi_v*Nvd+C_Phi*B_Phi_h], C_Delta, Csat);
D_Xi = [D_Phi_v*Nvh + D_Phi_h , zeros(4*n,2) ;...
        zeros(4,n) , D_Delta ;...
        Dsat1*Nuh , Dsat2 , zeros(4,1)];

eig_min = eigs(P_hat,1,0);
P_hat = eig_min \ P_hat; % normalize P_hat for numerical stability

cvx_begin sdp quiet
    cvx_solver Mosek
    variable P(n_Xi,n_Xi) symmetric
    variable muQ(n) nonnegative
    variable muH(n) nonnegative
    variable rho_v(n) nonnegative
    variable rho_h(n) nonnegative
    variable epsilon nonnegative
    variable eta_0(n) nonnegative

    variable sat_sec nonnegative
    variable sat_slop nonnegative
    variable rho_sat nonnegative
    variable rho_delta(n_Delta) nonnegative
    variable M11(2,2) semidefinite

    Px = P(1:n_G , 1:n_G);
    
    Qsec = [zeros(n,n) , diag(muQ) ; diag(muQ) , zeros(n,n)];
                 
    Qslop = [zeros(n) , diag(eta_0) ; diag(eta_0) , zeros(n)];

    M = blkdiag(delta_norm*M11, -1/delta_norm*M11);
    Qsat = blkdiag([0 , sat_sec ; sat_sec , 0] , [0 , sat_slop ; sat_slop , 0]);
              
    Ieps = blkdiag(epsilon*eye(n_G) , zeros(n_Xi-n_G+n+2));
    
    I0 = [eye(n_Xi) , zeros(n_Xi,n+2)];
    AB = [A_Xi , B_Xi];
    CD = [C_Xi , D_Xi];
    
    H = [-diag(rho_v./v_bound_square + 2*alpha_max.*alpha_min.*muH) , diag((alpha_max+alpha_min).*muH) ;...
          diag((alpha_max+alpha_min).*muH) , -diag(rho_h./h_bound_square + 2*muH)];
    
    minimize( trace(Px) )
    subject to
        I0' * P * I0 - AB' * P * AB >= Ieps + CD' * blkdiag(Qsec, Qslop, M, Qsat) * CD;
        P >= blkdiag((1+sum(rho_v+rho_h)+sum(rho_delta)+sat_bound_square*rho_sat)*P_hat , H/eig_min , ...
                      -diag(rho_delta./delta_bound_square)/eig_min, -rho_sat/eig_min);
cvx_end

Px = eig_min * Px;
if isnan(cvx_optval)
    if min(eig(I0' * P * I0 - AB' * P * AB - Ieps - CD' * blkdiag(Qsec, Qslop, M, Qsat) * CD)) >= 0 && ...
       min(eig(P - blkdiag((1+sum(rho_v+rho_h)+sum(rho_delta)+sat_bound_square*rho_sat)*P_hat , H/eig_min , ...
               -diag(rho_delta./delta_bound_square)/eig_min, -rho_sat/eig_min))) >= 0

    Px_trace = trace(Px);
else
    Px_trace = inf;
    end
else
     Px_trace = eig_min * cvx_optval;
end