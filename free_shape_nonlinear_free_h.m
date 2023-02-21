function out = ...
    free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, psi_dot_init)
% now solving the nonlinear shape equation using shooting method and
% fmincon
    min_psi_dot = -Inf;
    max_psi_dot = -1e-3;
    lambda = sqrt(kappa/Sigma);
    const(1) = lambda;
    const(2) = d;
    
    r_phi = R*sin(phi);
    psi_init = phi;
%     p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);
%     z = 1-cos(phi);
%     p_r_init = sqrt(z*(2-z))/(1-z)*(1/R+2*z*R/lambda^2-R*psi_dot_init^2);

    min_constants = [lambda, R, phi, d];

%     psi_end = pi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve using fmincon
%     options = optimoptions('fminunc', 'MaxFunctionEvaluations',1000)
    options = optimoptions('fmincon', 'algorithm','sqp', 'MaxFunctionEvaluations',1000);

%     outputs = fminsearch(@(inp) get_shape(inp, min_constants), [psi_dot_init]);
%     outputs = fminunc(@(inp) get_shape(inp, min_constants), psi_dot_init,...
%         options);
%     while psi_end>pi/2
    counter = 0;
    while true
        outputs = fmincon(@(inp) get_shape(inp, min_constants, 2), psi_dot_init,...
            [],[],[],[],min_psi_dot,max_psi_dot,[],options);
    
        psi_dot_init_final = outputs(1);
    
        [~,out] = get_shape(psi_dot_init_final, min_constants, 2);
        psi_end = out.ye(1);
        r_end = out.ye(2);
    
        if psi_end<0.001
            break
        end
    
        % if we don't have zero curvature at end, bracket
        % only do this the first time.
        if counter==0
            while psi_end>0.001||r_end>d/2
                counter = counter+1;
                if counter>100
                    break
                end
                if mod(counter,2)==0
                    psi_dot_init_new = psi_dot_init_final/1.2^counter;
                else
                    psi_dot_init_new = psi_dot_init_final*1.2^counter;
                end
                [diffs,out] = get_shape(psi_dot_init_new, min_constants, 2);
                psi_end = out.ye(1);
                r_end = out.ye(2);
                if psi_dot_init_new>psi_dot_init_final
                    min_psi_dot = psi_dot_init_final;
                    max_psi_dot = psi_dot_init_new;
                    max_correct = 1;
                else
                    min_psi_dot = psi_dot_init_new;
                    max_psi_dot = psi_dot_init_final;
                    max_correct = 0;
                end
                psi_dot_init_final = psi_dot_init_new;
                diffs;
            end
        else
            % otherwise, bracket the value
            if max_correct
                min_psi_dot = psi_dot_init;
            else
                max_psi_dot = psi_dot_init;
            end
%             if psi_end>0.001
%                 % bounds need to be more negative, upper bound is wrong!
%                 max_psi_dot = psi_dot_init_final;
%             elseif psi_end<=0.99*pi&&r_end>d/2
%                 % bounds need to be more negative, upper bound is wrong!
%                 max_psi_dot = psi_dot_init_final;
%             elseif psi_end<=0.99*pi&&r_end<d/2
%                 % bounds need to be less negative, upper bound is wrong!
%                 min_psi_dot = psi_dot_init_final;
%             end
        end
        
        if counter>100
            break
        end

        
        % bisect
        psi_dot_init = (min_psi_dot+max_psi_dot)/2;
%         counter = counter+1;
        max_psi_dot;
        min_psi_dot;

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve using bracketing method
    % start at initial guess, and hunt either up or down, where we want to
    % ultimately bracket the answer and then finally run fmincon with those
    % boundaries
    
%     psi_dot = psi_dot_init;
%     psi_dot_lower = psi_dot/2;
%     psi_dot_upper = psi_dot*2;
%     [diffs,out] = get_shape(psi_dot, min_constants, 4);
%     psi_end = out.ye(1);
%     r_end = out.ye(2);
%     % start by bracketing the output
%     if psi_end>0.99*pi
%         % if psi is too large, decrease until it's 0
%         while psi_end>0.99*pi
%             psi_dot_upper = psi_dot;
%             psi_dot = psi_dot/2;
%             [~,out] = get_shape(psi_dot, min_constants, 4);
%             psi_end = out.ye(1);
%             r_end = out.ye(2);
%         end
%         psi_dot_lower = psi_dot;
%     elseif psi_end<=0.99*pi&&r_end>d/2
%         % if psi pushes r past d/2, bracket upwards
%         while psi_end<=0.99*pi&&r_end>d/2
%             psi_dot_lower = psi_dot;
%             psi_dot = psi_dot*2;
%             [~,out] = get_shape(psi_dot, min_constants, 4);
%             psi_end = out.ye(1);
%             r_end = out.ye(2);
%         end
%         psi_dot_upper = psi_dot;
%     elseif psi_end<=0.99*pi&&r_end<d/2
%         while psi_end<=0.99*pi&&r_end>d/2
%             psi_dot_lower = psi_dot;
%             psi_dot = psi_dot*2;
%             [~,out] = get_shape(psi_dot, min_constants, 4);
%             psi_end = out.ye(1);
%             r_end = out.ye(2);
%         end
%         psi_dot_upper = psi_dot;
%     end
%     % now we have our brackets, we want to bisect 
%     while diffs>1e-10||psi_end>0.001*pi
%         if psi_end>0.001*pi
%             % if psi is too large, decrease until it's 0
%             psi_dot_upper = psi_dot;
%             psi_dot = psi_dot/2;
%             [~,out] = get_shape(psi_dot, min_constants, 4);
%             psi_end = out.ye(1);
%             r_end = out.ye(2);
%             diffs = (r_end-d/2)^2;
%             psi_dot_lower = psi_dot;
%         elseif psi_end<=0.99*pi&&r_end>d/2
%             % if psi pushes r past d/2, bracket upwards
%             while psi_end<=0.99*pi&&r_end>d/2
%                 psi_dot_lower = psi_dot;
%                 psi_dot = psi_dot*2;
%                 [~,out] = get_shape(psi_dot, min_constants, 4);
%                 psi_end = out.ye(1);
%                 r_end = out.ye(2);
%             end
%             psi_dot_upper = psi_dot;
%         elseif psi_end<=0.99*pi&&r_end<d/2
%             while psi_end<=0.99*pi&&r_end>d/2
%                 psi_dot_lower = psi_dot;
%                 psi_dot = psi_dot*2;
%                 [~,out] = get_shape(psi_dot, min_constants, 4);
%                 psi_end = out.ye(1);
%                 r_end = out.ye(2);
%             end
%             psi_dot_upper = psi_dot;
%         end
%     end






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init_final, p_r_init_final, 0,0,0];
%     event_func = @(s,y) myEventFcn(s,y,const);
% %     options = odeset('Events', @myEventFcn);
%     options = odeset('Events', event_func);
%     out = ode45(@(s, y) hamilton(s, y,const),[0,100],init_vals,...
%         options);
% 
% %         psi_end = out.ye(1);
% %         psi_dot_init = psi_dot_init*(rand-0.5)*5+rand/2
% % 
% %     end


end

function derivs = hamilton(s, y, const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);
    E_bend = y(7);
    A = y(8);

    lambda = const(1);
    d = const(2);

    f(1) = p_psi/(2*r)-sin(psi)/r;
    f(2) = cos(psi);
    f(3) = sin(psi);
    f(4) = cos(psi)*(p_psi/r-p_h)...
        +sin(psi)*p_r;
    f(5) = p_psi/r*(p_psi/(4*r)-sin(psi)/r)...
        +2/lambda^2;
    f(6) = 0;
    f(7) = p_psi^2/(4*r);
    f(8) = r;

    derivs = f';
end

function [value,isterminal,direction] = myEventFcn(s,y,const)
    psi = y(1);
    r = y(2);
    h = y(3);
    p_psi = y(4);
    p_r = y(5);
    p_h = y(6);

    lambda = const(1);
    d = const(2);

%     value = y(2)-d/2;
%     value = y(1);
%     isterminal = 1;
%     direction = 0;
% 
%     value = [y(1);y(1)-pi;y(2)-d/2];
%     isterminal = [1;1;1];
%     direction = [0;0;0];

    value = [y(1);y(1)-pi];
    isterminal = [1;1];
    direction = [0;0];

%     value = [y(1);y(2)-d/2];
%     isterminal = [1;1];
%     direction = [0;0];
end

function [diffs,out] = get_shape(psi_dot_init, constants, condition)

    lambda = constants(1);
    R = constants(2);
    phi = constants(3);
    d = constants(4);

    const(1) = lambda;
    const(2) = d;
    
    r_phi = R*sin(phi);
    psi_init = phi;
    p_psi_init = 2*r_phi*(psi_dot_init+sin(psi_init)/r_phi);
    p_r_init = -(p_psi_init^2/(4*r_phi)-p_psi_init*sin(psi_init)/r_phi...
        -2*r_phi/lambda^2)/cos(psi_init);

    init_vals = [psi_init, r_phi, R*sin(phi), p_psi_init, p_r_init, 0,0,0];
    event_func = @(s,y) myEventFcn(s,y,const);
%     options = odeset('Events', @myEventFcn);
    options = odeset('Events', event_func);
    out = ode45(@(s, y) hamilton(s, y,const),[0,100],init_vals,...
        options);
    
    if condition==1
        psi_end = out.y(1,end);
        diffs = psi_end^2;
    elseif condition==2
        r_end = out.y(2,end);
        diffs = (r_end-d/2)^2;
    elseif condition==3
        r_end = out.y(2,end);
        diffs = r_end-d/2;
    elseif condition==4   
        psi_end = out.y(1,end);
        r_end = out.y(2,end);
        diffs = (r_end-d/2)^2+50*psi_end^2;
    end

%     if psi_end>pi/2
%         diffs = diffs*1000;
%     end
end