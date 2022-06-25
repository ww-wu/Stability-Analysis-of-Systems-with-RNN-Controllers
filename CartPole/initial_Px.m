function Px = initial_Px(sys,Whx,Whh)

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

A = [A_G , zeros(n_G,n) ; zeros(n,n_G), zeros(n)];
B = [B_G*Nuh ; eye(n)];

K = (eye(n)-Nvh)\[Nvy*C_G , Nvd];

cvx_begin sdp quiet
    cvx_solver Mosek
    variable P(n_G+n,n_G+n) symmetric
    
    Px = P(1:n_G , 1:n_G);

    maximize( 0 )
    subject to
        P - (A+B*K)'*P*(A+B*K) >= 0;
        P >= eye(n_G+n);
cvx_end