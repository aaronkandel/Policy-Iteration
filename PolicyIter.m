function optimalPolicy = PolicyIter(stateTrans,stateCost,nInputs,nStates,discountFactor)
% Penn State ME597.001: Optimal Control of Energy Systems

% Begin by choosing an arbitrary control policy
currentPolicy = ones(nStates,1); % 
previousPolicy = currentPolicy;

% Next, choose an arbitary initial guesstimate for the value function
valueFunction = zeros(nStates,1);
policyError = 5;   % Difference in control actions between two consecutively optimized control policies

while policyError > 0.00005     % Need policy error to eventually hit zero 
    valueFunction = IterPolEval(stateTrans,stateCost,...
    nInputs,nStates,currentPolicy,discountFactor,0.000001,valueFunction);
    
    bestValueFunction = valueFunction; 
    for currentStateIndex = 1:nStates
        for controlIndex = 1:nInputs
            currentValueFunction = 0;
            nextStateIndex = stateTrans(currentStateIndex,controlIndex);
            currentValueFunction = currentValueFunction + ...
                    discountFactor*1*(valueFunction(nextStateIndex) + stateCost(currentStateIndex,controlIndex));
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

