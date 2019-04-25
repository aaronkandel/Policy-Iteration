function optimalPolicy = PolicyIterProb(stateTrans,stateCost,nInputs,nStates,discountFactor)
% Penn State ME597.001: Optimal Control of Energy Systems

% Begin by choosing an arbitrary control policy
currentPolicy = ones(nStates,1); % 
previousPolicy = currentPolicy;

% Next, choose an arbitary initial guesstimate for the value function

valueFunction = zeros(nStates,1);

policyError = 5;   % Difference in control actions between two consecutively optimized control policies

while policyError > 0.00005     % Need policy error to eventually hit zero 
    valueFunction = IterPolEvalProb(stateTrans,stateCost,...
    nInputs,nStates,currentPolicy,discountFactor,0.000001,valueFunction);
    
    
    bestValueFunction = valueFunction; 
    for currentStateIndex = 1:nStates
        for controlIndex = 1:nInputs
            currentValueFunction = 0;
            for nextStateIndex = 1:nStates          
                currentValueFunction = currentValueFunction + ...
                    discountFactor*markovTransitionTables(currentStateIndex,nextStateIndex,controlIndex)*(valueFunction(nextStateIndex) + ...
                    markovTransitionCosts(currentStateIndex,nextStateIndex,controlIndex));
            end
            if currentValueFunction < bestValueFunction(currentStateIndex)
                bestValueFunction(currentStateIndex) = currentValueFunction;
                currentPolicy(currentStateIndex) = controlIndex;
            end
        end
    end
    policyError = sum(abs(currentPolicy-previousPolicy));
    previousPolicy = currentPolicy;

end 

optimalPolicy = currentPolicy;

end

