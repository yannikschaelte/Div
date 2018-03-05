function [ list ] = predict_outliers(data,simfun,n_runs,max_fun_evals)
% do n_runs (quick) multistart optimizations with budgets of max_fun_evals
% function evaluations, and then use the obtained parameters for a
% comparison with the data. extract parameters which contribute most to the
% negative log-posterior


