rng(298)
tic

%specify the number of processes and goods

global n m A B r_total r_count;
n = 100;        %process
m = 15;         %good

%%A and B are positive
A = abs(0.1 * randn(n,m) + 1)+1e-5;                  %consumption function
B = abs(0.1 * randn(n,m) + 1)+1e-5;      
      
%production function
r0 = ones(n,1);                                      %initial guess of intensity function
q0 = 10 + randn(m,1);                                %initial guess of q 
q0 = q0/sqrt(sum(q0.^2));                            %normalize q

%set up for simulation
t_max = 100;                        %max timestep
q_total = zeros(m,t_max + 1);       %record q
p_total = zeros(m,t_max + 1);       %record p
excess_total = zeros(m,t_max + 1);  %record excess demand
r_total = zeros(n,t_max + 2);       %record intensity
r_total(:,1) = r0;              
r_count = 1;                    %t is not global, can't be achieved inside the function, so use r_count increase together
                                %with t to keep reference

fun = @Minimization;            %objective function to minimize

nonlicon = @nlcon;              %nonlinear constraint
options = optimoptions('fmincon',"EnableFeasibilityMode",true,'SpecifyObjectiveGradient',true, ...
    'MaxIterations', 3e3, 'MaxFunctionEvaluations',3e4, "SubproblemAlgorithm","cg",'Display','iter');   %options for fmincon

loop = zeros(1,t_max+1);          %records the # of loops we need to find a minimum                         
for t = 0:t_max 
    if ((t+1) == 1)
        q_equi = fmincon(fun, q0, [], [], [], [], [], [], nonlicon, options);  %minimize the objective function with equilibiurm price in the previous day
    else
        q_equi = fmincon(fun, q_total(:,t), [], [], [], [], [], [], nonlicon, options);  %minimize the objective function   
    end
    % if fun(q_equi) > 0.01 (which should be zero theoretically, but we have some tolerance in practice),
    % then q_equi is not a global minimizer, then we rerun the minimization
    % with a different starting point
    while fun(q_equi) > 0.01
        q0 = 10 + randn(m,1);                      %initial guess of q 
        q0 = q0/sqrt (sum(q0.^2));                 %normalize q
        q_equi = fmincon(fun, q0, [], [], [], [], [], [], nonlicon, options);  %minimize the objective function
        loop(1,t+1) = loop(1,t+1) + 1;
    end
    p_equi = q_equi.^2;             %equilibrium price
    %record in correspond matrices
    q_total(:,t+1) = q_equi; 
    p_total(:,t+1) = p_equi;
    excess = ExcessDemand(q_equi);
    
    for count = 1:m         %eliminate small excess demand
        if (excess(count) < 1e-2 && excess(count) > 0) || (excess(count) > -1e-2 && excess(count) < 0)
            excess_total(count,t+1) = 0;
        else
            excess_total(count,t+1) = excess(count);
        end
    end

    %compute r(t+1) based on r(t)
    r_total(:,(t+1)+1) = (r_total(:,t+1) .* ((B*(q_equi.^2)) ./ (A*(q_equi.^2)))) ; 
    for count = 1:n         %eliminate small intensity
        if (r_total(count,(t+1)+1)) < 1e-2 && (r_total(count,t+1)) > 0 || r_total(count,(t+1)+1) < 0
            r_total(count,(t+1)+1) = 0;
        end
    end

    r_count = r_count + 1;      %alternative referencing t 
end

toc

%Defined Functions
function [E, grad_E] = ExcessDemand(q)
    global m A B r_total r_count;
    %excess demand function
    E = ( transpose(r_total(:,r_count)) * (repmat((B*(q.^2)) ./ (A*(q.^2)), 1, m) .* A - B) ).';
    %gradient of excess demand function
    grad_E = 2.*q.*((transpose(r_total(:,r_count)) * ((A.*( ((B.*repmat((A*(q.^2)), 1, m)) -(A.*repmat((B*(q.^2)), 1, m))) ./ (repmat((A*(q.^2)).^2, 1, m)) )))))';
end

%construct minimization function
function [phi, grad] = Minimization(q)
    global m A B r_total r_count;
    E = ( transpose(r_total(:,r_count)) * ( repmat((B*(q.^2)) ./ (A*(q.^2)), 1, m) .* A - B ) ).'; %excess demand function
    grad_E = 2.*q.*((transpose(r_total(:,r_count)) * ((A.*( ((B.*repmat((A*(q.^2)), 1, m)) -(A.*repmat((B*(q.^2)), 1, m))) ./ (repmat((A*(q.^2)).^2, 1, m)) )))))';

    phi = 0;
    grad = zeros(m,1);
    for e_j = 1:m                       %constructing minimization function 
        if E(e_j) > 0
            phi = phi + E(e_j);
        end
    end

    %compute gradient of excess dmeand function
    if nargout > 1 % gradient required
        for k = 1:m
            if (E(k) > 0) 
                grad(k) = grad(k) + (E(k) * (grad_E(k))') ;
            end
        end
    end
end

function [c,ceq] = nlcon(q)
    c = [];
    ceq = sum(q.^2)-1;
end





%debug area
% function [E, Sum_matrix] = teste(q)
%     global n m A B r_total r_count;
%     E = zeros(m,1);
%     for k = 1:m
%         E_n=0;
%         Sum_matrix = zeros(n,1);
%         for count_n = 1:n
%             Bn=0;
%             An=0;
%             for count_m = 1:m
%                 Bn = Bn + B(count_n,count_m) * q(count_m)^2;
%                 An = An + A(count_n,count_m) * q(count_m)^2;
%             end
%             Sum_matrix(count_n) = (Bn/An*A(count_n,k)-B(count_n,k));
%             E_n = E_n + (Bn/An*A(count_n,k)-B(count_n,k));
%         end
%         E(k) = E_n;
%     end
% end
