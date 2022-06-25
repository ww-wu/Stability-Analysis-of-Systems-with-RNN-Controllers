clear
load('lateral_vehicle.mat')

P_hat = initial_Px(sys,Whx,Whh);
while true
    [Px,optval] = ROA_sat(sys, Whx, Whh, P_hat);
    if isinf(optval)
        P_hat = P_hat * 10;
    else
        break
    end
end
Pxx = Px;
P_min = zeros(size(Px));
P_max = P_hat;
P_hat = (Px + P_min) / 2;

T = 40;
trace_hat = zeros(1,T);
trace_Px = zeros(1,T);
trace_min = zeros(1,T);
trace_max = zeros(1,T);


for k = 1:T
    [Px,optval] = ROA_sat(sys, Whx, Whh, P_hat);
    disp(optval)

    trace_hat(k) = trace(P_hat);
    trace_Px(k) = optval;
    trace_min(k) = trace(P_min);
    trace_max(k) = trace(P_max);

    if optval < trace(Pxx)
        Pxx = Px;
        P_max = P_hat;
    else
        if isinf(optval) || max(eig(P_hat - P_max)) < 0
            P_min = P_hat;
        else
            P_hat =(P_max + P_min)/2;
            continue
        end
    end

    if eigs(P_max-P_min,1) < eigs(P_min,1,0)
        break
    end
    P_hat = (Pxx + P_min)/2;
end
