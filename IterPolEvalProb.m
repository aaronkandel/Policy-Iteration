function valueFunction = IterPolEvalProb(stateTrans,stateCost,...
    nInputs,nStates,controlPolicy,discountFactor,tolerance,initialEstimate)
% Penn State ME597.001: Optimal Control of Energy Systems

% Begin by setting initial estimate for value function
valueFunction = initialEstimate; 
newValueFunction = valueFunction;

% Next, construct the overall Markov transition probability table, and
% overall Markov transition cost table

overallMarkovTransitionTable = zeros(nStates,nStates);
overallMarkovTransitionCostTable = zeros(nStates,nStates);
for stateIndex = 1:nStates
    controlIndex = controlPolicy(stateIndex);
    overallMarkovTransitionTable(stateIndex,:) = squeeze(markovTransitionTables(stateIndex,:,controlIndex));%squeeze(markovTransitionTables(stateIndex,:,controlIndex));
    overallMarkovTransitionCostTable(stateIndex,:) = squeeze(markovTransitionCosts(stateIndex,:,controlIndex));%squeeze(markovTransitionCosts(stateIndex,:,controlIndex));
end

% Next, perform value iteration
valueError = 5*tolerance; % 10

while valueError > tolerance
    for k = 1:nStates
        newValueFunction(k) = 0;
        for i = 1:nStates
            newValueFunction(k) = newValueFunction(k) + discountFactor*...
                overallMarkovTransitionTable(k,i)*(valueFunction(i)+...
                overallMarkovTransitionCostTable(k,i));
        end
    end
    valueError = (sqrt(sum((newValueFunction-valueFunction).^2))/sqrt(sum(valueFunction.^2)))*100;
    valueFunction = newValueFunction;
end
end

