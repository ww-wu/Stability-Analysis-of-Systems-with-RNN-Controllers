function [v_max,v_min,h_max,h_min,alpha_max,alpha_min,sigma_max,sigma_min] = interval_bound_propagation(C_G,Whx,Whh,P)

L = length(Whh);
n_G = size(C_G,2);
xs = zeros(n_G,1); % Steady State

vs = cell(L,1); % Pre-activation Vector associated to Steady State
v_max = cell(L,1);
v_min = cell(L,1);

hs = cell(L,1); % Hidden State associated to Steady State
h_max = cell(L,1);
h_min = cell(L,1);

alpha_max = cell(L,1); % Sector Upper Bound
alpha_min = cell(L,1); % Sector Lower Bound
sigma_max = cell(L,1); % Slope Upper Bound
sigma_min = cell(L,1); % Slope Lower Bound

for e = 1:L
    if e == 1
        dv = sqrt(diag((Whx{e}*C_G) / P * (Whx{e}*C_G)'));
        vs{e} = Whx{e} * C_G * xs;
    else
        vs{e} = Whx{e} * hs{e-1};
        dv = abs(Whx{e}) * h_max{e-1};
    end
    Wabs = abs(Whh{e});
    v_max1 = dv + Wabs*ones(size(dv)); % +b{e} if nonzero
    v_min1 = -dv - Wabs*ones(size(dv));  % +b{e} if nonzero
    while true
        v_max2 = v_max1;
        v_max1 = dv + Wabs*tanh(v_max2); % +b{e} if nonzero

        v_min2 = v_min1;
        v_min1 = -dv - Wabs*tanh(v_min2); % +b{e} if nonzero
        if mean(v_max2 - v_max1) < 1e-3 && mean(v_min1 - v_min2) < 1e-3
            break
        end
    end
    v_max{e} = v_max1;
    v_min{e} = v_min1;

    hs{e} = tanh(vs{e});
    h_max{e} = tanh(v_max{e});
    h_min{e} = tanh(v_min{e});

    % Sector Bound
    a_max = (h_max{e}-hs{e}) ./ (v_max{e}-vs{e});
    a_min = (h_min{e}-hs{e}) ./ (v_min{e}-vs{e});

    % Slope Bound
    d_max = 1 - h_max{e}.^2; % derivative of tanh at v_max
    d_min = 1 - h_min{e}.^2; % derivative of tanh at v_min

    % Case 1: The input interval lies in the negative part, where tanh is convex
    flag1 = v_max{e} <= 0;
    alpha_max{e}(flag1) = a_max(flag1);
    alpha_min{e}(flag1) = a_min(flag1);
    sigma_max{e}(flag1) = d_max(flag1);
    sigma_min{e}(flag1) = d_min(flag1);

    % Case 2: The input interval lies in the positive part, where tanh is concave
    flag2 = v_min{e} >= 0;
    alpha_max{e}(flag2) = a_min(flag2);
    alpha_min{e}(flag2) = a_max(flag2);
    sigma_max{e}(flag2) = d_min(flag2);
    sigma_min{e}(flag2) = d_max(flag2);

    % Case 3: The origin is in the input interval
    flag3 = ~(flag1 | flag2);
    alpha_min{e}(flag3) = min(a_min(flag3), a_max(flag3));
    alpha_max{e}(flag3 & (vs{e}==0)) = 1;

    % for input inverals with nonzero vs, determine the sector upper bound
    % by using the bisection method to find the corresponding tangent point
    flag4 = flag3 & (vs{e}~=0);
    if any(flag4)
        % If vs > 0, then the tangent point is in [v_min,0]
        flag4_pos = vs{e}(flag4) > 0;
        v_sec_min = flag4_pos .* v_min{e}(flag4);
        % If vs < 0, then the tangent point is in [0,v_max]
        flag4_neg = vs{e}(flag4) < 0;
        v_sec_max = flag4_neg .* v_max{e}(flag4);
        while any(v_sec_max - v_sec_min > 1e-4)
            v_sec_mid = (v_sec_max + v_sec_min) / 2;
            h_sec_mid = tanh(v_sec_mid);
            % The slope of the secant intersecting the points (v_sec_mid,h_sec_mid) and (vs,hs)
            slop_sec = (h_sec_mid - hs{e}(flag4)) ./ (v_sec_mid - vs{e}(flag4));
            dh_sec_mid = 1 - h_sec_mid.^2; % derivative of tanh at v_sec_mid

            % v_sec_mid is further away the origin than the tangent point
            flag1 = slop_sec > dh_sec_mid;
            v_sec_min(flag1 & flag4_pos) = v_sec_mid(flag1 & flag4_pos);
            v_sec_max(flag1 & flag4_neg) = v_sec_mid(flag1 & flag4_neg);

            % v_sec_mid is closer to the origin than the tangent point
            flag2 = dh_sec_mid > slop_sec;
            v_sec_max(flag2 & flag4_pos) = v_sec_mid(flag2 & flag4_pos);
            v_sec_min(flag2 & flag4_neg) = v_sec_mid(flag2 & flag4_neg);
        end
        alpha_max{e}(flag4) = slop_sec;
    end

    sigma_max{e}(flag3) = 1;
    sigma_min{e}(flag3) = min(d_min(flag3), d_max(flag3));
end