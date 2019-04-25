function valueFunction = IterPolEval(stateTrans,stateCost,...
    nInputs,nStates,controlPolicy,discountFactor,tolerance,initialEstimate)
% Penn State ME597.001: Optimal Control of Energy Systems

% Begin by setting initial estimate for value function
valueFunction = initialEstimate; 
newValueFunction = valueFunction;

% Next, construct the overall Markov transition probability table, and
% overall Markov transition cost table
ovTransTable = zeros(nStates,1);
ovTransCostTable = zeros(nStates,1);

for stateIndex = 1:nStates
    controlIndex = controlPolicy(stateIndex);
    ovTransTable(stateIndex) = squeeze(stateTrans(stateIndex,controlIndex));
    ovTransCostTable(stateIndex) = squeeze(stateCost(stateIndex,controlIndex));
end

% Next, perform value iteration
valueError = 5*tolerance; % 10

while valueError > tolerance
    for k = 1:nStates
        newValueFunction(k) = 0;
        newValueFunction(k) = newValueFunction(k) + discountFactor*...
            1*(valueFunction(k)+... % ovTransTable(k)
            ovTransCostTable(k));%stateCost(k,u));
    end
    valueError = (sqrt(sum((newValueFunction-valueFunction).^2))/sqrt(sum(valueFunction.^2)))*100;
    valueFunction = newValueFunction;
end
end

