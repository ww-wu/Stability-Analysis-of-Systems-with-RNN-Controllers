function [Px, Px_trace] = ROA_full(sys,Whx,Whh,P_hat)

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

B_Psi_v = [-max(sigma_max)*eye(n) ; min(sigma_min)*eye(n)];
B_Psi_h = [eye(n) ; -eye(n)];
C_Psi = [zeros(2*n) ; zeros(2*n) ; eye(2*n)];
D_Psi_v = [diag(alpha_max) ; -diag(alpha_min);...
           max(sigma_max)*eye(n) ; -min(sigma_min)*eye(n) ; zeros(2*n,n)];
D_Psi_h = [-eye(n) ; eye(n) ; -eye(n) ; eye(n) ; zeros(2*n,n)];

n_Xi = n_G + 2*n;
A_Xi = [A_G , zeros(n_G,2*n) ; Nvy*C_G , zeros(n,n) , Nvd ; zeros(n,n_Xi)];
B_Xi = [B_G*Nuh ; Nvh ; eye(n)];
C_Xi = [D_Psi_v*Nvy*C_G , C_Psi*B_Psi_v , D_Psi_v*Nvd+C_Psi*B_Psi_h];
D_Xi = D_Psi_v*Nvh + D_Psi_h;

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
    variables eta_0(n,n) eta_n(n,n) eta_p(n,n)
    
    Px = P(1:n_G , 1:n_G);
    
    Qsec = [zeros(n,n) , diag(muQ) ; diag(muQ) , zeros(n,n)];
    
             
    Qslop = [zeros(n) , eta_0' , zeros(n), eta_p';
                  eta_0 , zeros(n) , eta_n, zeros(n);
                  zeros(n) , eta_n' , zeros(n,2*n);
                  eta_p, zeros(n,3*n)];
              
    Ieps = blkdiag(epsilon*eye(n_G) , zeros(n_Xi-n_G+n));
    
    I0 = [eye(n_Xi) , zeros(n_Xi,n)];
    AB = [A_Xi , B_Xi];
    CD = [C_Xi , D_Xi];
    
    H = [-diag(rho_v./v_bound_square + 2*alpha_max.*alpha_min.*muH) , diag((alpha_max+alpha_min).*muH) ;...
          diag((alpha_max+alpha_min).*muH) , -diag(rho_h./h_bound_square + 2*muH)];
    
    minimize( trace(Px) )
    subject to
        I0' * P * I0 - AB' * P * AB >= Ieps + CD' * blkdiag(Qsec, Qslop) * CD;
        P >= [(1+sum(rho_v+rho_h))*P_hat, zeros(n_G,2*n) ; zeros(2*n,n_G) , H/eig_min];
        for k = 1:n
            2*eta_0(k,k) >= sum(abs(eta_0(k,:))) + sum(abs(eta_n(k,:))) + sum(abs(eta_p(k,:)));
            2*eta_0(k,k) >= sum(abs(eta_0(:,k))) + sum(abs(eta_n(:,k))) + sum(abs(eta_p(:,k)));
        end
cvx_end

Px = eig_min * Px;
if isnan(cvx_optval)
    feasible = true;
    for k = 1:n
        if 2*eta_0(k,k) < sum(abs(eta_0(k,:))) + sum(abs(eta_n(k,:))) + sum(abs(eta_p(k,:)))
            feasible = false;
        end
        if 2*eta_0(k,k) < sum(abs(eta_0(:,k))) + sum(abs(eta_n(:,k))) + sum(abs(eta_p(:,k)))
            feasible = false;
        end
    end
    if min(eig(I0' * P * I0 - AB' * P * AB - Ieps - CD' * blkdiag(Qsec_diag,Qslop_diag) * CD)) < 0
        feasible = false;
    end
    if min(eig(P - [(1+sum(rho_v+rho_h))*P_hat, zeros(n_G,2*n) ; zeros(2*n,n_G) , H/eig_min])) < 0
        feasible = false;
    end

    if feasible
        Px_trace = trace(Px);
    else
        Px_trace = inf;
    end
else
    Px_trace = eig_min * cvx_optval;
end