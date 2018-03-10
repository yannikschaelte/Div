classdef NoodleResults < handle

    properties ( GetAccess = 'public', SetAccess = 'private' )
        x0;
        final_x;
        final_fval;
        final_grad;
        final_hess;
        final_gradnorm;
        iter_count;
        feval_count;
        cpu_time;
        exitflag;
    end
    
    methods
        
        function results = NoodleResults(problem)
            % Constructor
            if isa(problem,'noodles.NoodleProblem')
                results.x0              = problem.x0;
                results.final_x         = problem.state.x;
                results.final_fval      = problem.state.fval;
                results.final_grad      = problem.state.grad;
                results.final_hess      = problem.state.hess;
                results.final_gradnorm  = problem.state.gradnorm;
                results.iter_count      = problem.state.iter_count;
                results.feval_count     = problem.state.feval_count;
                results.cpu_time        = cputime - problem.start_time;
                results.exitflag        = problem.exitflag;
            else
                error('The NoodleResults constructor needs a NoodleProblem for initialization.');
            end
        end
        
    end
    
end

