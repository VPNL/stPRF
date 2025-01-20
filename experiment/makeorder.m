function order = makeorder(numconds,trialspercond)
% Generates counterbalanced order of conditions.
% Written by KGS Lab
% Edited by AS 8/2014

% check that numconds is a factor of trialspercond
numtrials = numconds*trialspercond;
if (numtrials/(numconds*numconds)) ~= round((numtrials/(numconds*numconds)))
    error('Total number of trials be a multiple of numconds.')
end

% set up goal state
order = zeros(numtrials,1);
for cond = 1:numconds
    for t = 1:trialspercond
        order(((cond-1)*trialspercond)+t) = cond;
    end
end
order = Shuffle(order);
goal = ones(numconds,numconds)*(numtrials/(numconds*numconds));

% minimize difference between the goal and the current history
while 1
    % get the energy for the current design
    history = gethistory(order);
    oldenergy = sum(sum(abs(history-goal)));
    % make a random swap in a copy of the design
    new = order;
    a = randi(numtrials);
    b = randi(numtrials);
    swap = new(a);
    new(a) = new(b);
    new(b) = swap;
    % calculate and minimize the energy in the new design
    newenergy = sum(sum(abs(gethistory(new)-goal)));
    if newenergy > 13
        if newenergy < oldenergy
            order = new;
        end
    else
        if newenergy <= oldenergy
            order = new;
        end
    end
    if newenergy == 1
        break
    end
end

end