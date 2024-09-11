function [fMin, bestX, Convergence_curve] = SLPDBO(pop, M, c, d, dim, fobj)

    P_percent = 0.2;    % The population size of producers accounts for "P_percent" percent of the total population size
    pNum = round(pop * P_percent);    % The population size of the producers   初始化加莱维飞行

    lb = c .* ones(1, dim);    % Lower limit/bounds/ a vector
    ub = d .* ones(1, dim);    % Upper limit/bounds/ a vector

%     % Initialization using chaotic mapping
%     x = generateChaoticPopulation(pop, dim, lb, ub);
%     for i = 1:pop
%         fit(i) = fobj(x(i, :)',func_num );
%     end
%% 改进初始化阶段--加入sinusoidal映射
        sinusoidal=2.3;
        Sinusoidal=rand(pop,dim);
        for i=1:pop
            for j=2:dim
                Sinusoidal(i,j)=sinusoidal*Sinusoidal(i,j-1).^2*(sin(pi*Sinusoidal(i,j-1)));
            end
        end
        result = Sinusoidal;
        x=repmat(lb,pop,1)+result.*repmat((ub-lb),pop,1);
        for i=1:pop
            fit( i ) = fobj( x( i, : ) ) ; 
        end
    pFit = fit;
    pX = x;
    XX = pX;    
    [fMin, bestI] = min(fit);      % fMin denotes the global optimum fitness value
    bestX = x(bestI, :);           % bestX denotes the global optimum position corresponding to fMin

    % PSO parameters
    w_max = 0.9;
    w_min = 0.4;
    c1 = 1.5; % Cognitive parameter
    c2 = 1.5; % Social parameter
    v = zeros(pop, dim); % Velocity initialization

    % Start updating the solutions
    for t = 1:M
        display(['The iter is ', num2str(t)]);
        [fmax, B] = max(fit);
        worse = x(B, :);   
        r2 = rand(1);
        
        R = 1 - t / M; % Adaptive parameter
        
        % Adaptive inertia weight
        w = w_max - ((w_max - w_min) * t / M);
        
        % Adaptive mutation probability
        mutation_prob = 0.1 * (1 - t / M); % Decrease over time

        % Producers update
        for i = 1:pNum
            if r2 < 0.9
                r1 = rand(1);
                a = rand(1, 1);
                if a > 0.1
                    a = 1;
                else
                    a = -1;
                end
                x(i, :) = pX(i, :) + 0.3 * abs(pX(i, :) - worse) + a * 0.1 * (XX(i, :)); % Equation (1)
            else
                aaa = randperm(180, 1);
                if aaa == 0 || aaa == 90 || aaa == 180
                    x(i, :) = pX(i, :);
                end
                theta = aaa * pi / 180;
                x(i, :) = pX(i, :) + tan(theta) .* abs(pX(i, :) - XX(i, :)); % Equation (2)
            end
            x(i, :) = Bounds(x(i, :), lb, ub);
            fit(i) = fobj(x(i, :) );
        end

        [fMMin, bestII] = min(fit);      % fMin denotes the current optimum fitness value
        bestXX = x(bestII, :);           % bestXX denotes the current optimum position 

        Xnew1 = bestXX .* (1 - R); 
        Xnew2 = bestXX .* (1 + R);       %%% Equation (3)
        Xnew1 = Bounds(Xnew1, lb, ub);
        Xnew2 = Bounds(Xnew2, lb, ub);

        Xnew11 = bestX .* (1 - R); 
        Xnew22 = bestX .* (1 + R);       %%% Equation (5)
        Xnew11 = Bounds(Xnew11, lb, ub);
        Xnew22 = Bounds(Xnew22, lb, ub);

        % Scroungers update
        for i = (pNum + 1):12              % Equation (4)
            x(i, :) = bestXX + ((rand(1, dim)) .* (pX(i, :) - Xnew1) + (rand(1, dim)) .* (pX(i, :) - Xnew2));
            x(i, :) = Bounds(x(i, :), Xnew1, Xnew2);
            fit(i) = fobj(x(i, :));
        end

        for i = 13:19                      % Equation (6)
            x(i, :) = pX(i, :) + ((randn(1)) .* (pX(i, :) - Xnew11) + ((rand(1, dim)) .* (pX(i, :) - Xnew22)));
            x(i, :) = Bounds(x(i, :), lb, ub);
            fit(i) = fobj(x(i, :) );
        end

        for j = 20:pop                     % Equation (7)
            RL = 0.15*levy(pop,dim,1.5);
            %x(j, :) = bestX + randn(1, dim) .* ((abs((pX(j, :) - bestXX))) + (abs((pX(j, :) - bestX)))) / 2;
            x( j,: )=bestX+(randn(1,dim).*((abs(( pX(j,:  )-bestXX)))+(abs(( pX(j,:  )-bestX))))./2).*RL(i,:);  %加入莱维飞行
            x(j, :) = Bounds(x(j, :), lb, ub);
            fit(j) = fobj(x(j, :) );
        end

        % PSO velocity and position update
        for i = 1:pop
            if r2 < 0.9
                v(i, :) = w * v(i, :) + c1 * rand(1, dim) .* (pX(i, :) - x(i, :)) + c2 * rand(1, dim) .* (bestX - x(i, :));
                x(i, :) = x(i, :) + v(i, :);
            end

            % Apply adaptive mutation
            if rand < mutation_prob
                x(i, :) = x(i, :) + 0.1 * randn(1, dim);
            end

            
            x(i, :) = Bounds(x(i, :), lb, ub);
            fit(i) = fobj(x(i, :));
        end

        % Update the individual's best fitness value and the global best fitness value
        XX = pX;
        for i = 1:pop
            if fit(i) < pFit(i)
                pFit(i) = fit(i);
                pX(i, :) = x(i, :);
            end
            
            if pFit(i) < fMin
                fMin = pFit(i);
                bestX = pX(i, :);
            end
        end
        
        Convergence_curve(t) = fMin;
    end
end

% Application of simple limits/bounds
function s = Bounds(s, Lb, Ub)
    % Apply the lower bound vector
    temp = s;
    I = temp < Lb;
    temp(I) = Lb(I);
    
    % Apply the upper bound vector 
    J = temp > Ub;
    temp(J) = Ub(J);
    % Update this new move 
    s = temp;
end
%%
function [z] = levy(n,m,beta)
num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator
den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator
sigma_u = (num/den)^(1/beta);% Standard deviation
u = random('Normal',0,sigma_u,n,m);
v = random('Normal',0,1,n,m);
z =u./(abs(v).^(1/beta));
end
% % Initial population generation using chaotic mapping
% function chaoticPopulation = generateChaoticPopulation(pop, dim, lb, ub)
%     chaoticSequence = zeros(pop, dim);
%     x = rand(); % Initial value for chaotic sequence
%     for i = 1:pop
%         for j = 1:dim
%             x = 4 * x * (1 - x); % Logistic map
%             chaoticSequence(i, j) = lb(j) + (ub(j) - lb(j)) * x;
%         end
%     end
%     chaoticPopulation = chaoticSequence;
% end

