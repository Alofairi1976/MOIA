function [ MDS ] = GRBMDS_ILP(E )
% This function is to solve ILP for Minimum Dominating set problem
% Using gurobi optimization solver.
tic
n=size(E,1);

    clear model;
    model.A = sparse(E);
    model.obj = ones(n,1);
    model.rhs = ones(n,1);
    model.lb = zeros(n,1);
    model.ub = ones(n,1);
    model.sense = '>' ;
    model.vtype = 'B';
    model.modelsense = 'min';
 
    clear params;
    params.outputflag = 0;
    gurobi_write(model, 'sos.lp');

    result = gurobi(model, params);
    MDS=find(result.x); 
   
toc
end

