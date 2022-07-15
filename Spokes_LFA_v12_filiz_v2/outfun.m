function stop = outfun(x, optimValues, state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
stop = false;

global firstorderopt
global time

switch state
    case 'init'
        tic
    case 'iter'
          tmp = toc;
          time = [time; tmp];
          tic
          firstorderopt = [firstorderopt; optimValues.firstorderopt];
    case 'done'
        toc
    otherwise
end

end