function[result, f] = alg2g(problem)
% Solves a  Mixed-Integer Linear Optimization Problem.
% Input
%   problem - structure of the problem.
% Output
%   xopt - Solution vector.
%   fopt - Solution function value.
%   itcn - iteration count

%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64')
%addpath('C:\gurobi910\win64\matlab')

% LP problem
[meq, neq] = size(problem.Aeq);
[mineq, nineq] = size(problem.Aineq);
n=neq;

clear model;
model.A = sparse([problem.Aeq;problem.Aineq]);
model.obj = problem.f;
model.rhs = [problem.beq;problem.bineq];
sense = char(1,meq+mineq);
for i = 1:meq
    sense = '=';
end
for i = meq+1:meq+mineq
    sense = '<';
end
model.sense = sense;
model.ub = problem.ub;
model.lb = problem.lb;
model.modelsense = 'min';

result = gurobi(model);

Aineq = [problem.Aineq, zeros(mineq,1+3*n);
         zeros(mineq,nineq), problem.Aineq, zeros(mineq,1+2*n)];

bineq = [problem.bineq; problem.bineq];
beq = [result.x(1:n); zeros(n,1)];

alpha = 0.5*max(abs(problem.f));
beta = 5*(alpha)*sign(problem.f);

model.obj = [problem.f; zeros(n,1); alpha; -beta; beta];

% Find the basis
A = [problem.Aineq, eye(mineq);
    problem.Aeq, zeros(meq,mineq)];

base = [];

for i = 1:n
    if result.vbasis(i) == 0
        base = [base; i];
    end
end
baseorig = base;
for i = 1:mineq
    if result.cbasis(i) == 0
        base = [base; i+n];
    end
end

BR  = A(:, base);
nbase = setdiff(1:n+mineq,base);
NBR = A(:, nbase);

[~,nnb] = size(NBR);
D = zeros(n,nnb);

for i = 1:nnb
    aux = BR\NBR(:,i);
    for j = 1:length(baseorig)
        D(baseorig(j),i) = aux(j);
    end
    if nbase(i) < n
        D(nbase(i),i) = -1;
    end
end


% Generates the equivalent model:
%  Variables order: (y, x, lambda, p).

%for it = 1:nnb
    
    c = [problem.f;zeros(mineq,1)];
    
    y = BR'\c(base);
    rn = c(nbase) - NBR'*y;
    
    [~,it] = max(rn);

    d = D(problem.intcon,it);
    eyeN = eye(n);

    sd = length(d);
    
    Aeq = [problem.Aeq,            zeros(meq,n),           zeros(meq,1), zeros(meq,n),              zeros(meq,n);
            zeros(sd,n),            eyeN(problem.intcon,:), d,             zeros(sd,n),             zeros(sd,n);
            eyeN(problem.intcon,:), -eyeN(problem.intcon,:), zeros(sd,1), eyeN(problem.intcon,:),-eyeN(problem.intcon,:)];
    
    x = result.x(1:n);

    beq = [problem.beq; x(problem.intcon); zeros(sd,1)];
    
    model.lb = [problem.lb; problem.lb; 0; zeros(n,1);zeros(n,1)];
    model.ub = [problem.ub; problem.ub; inf; inf(n,1); inf(n,1)];
    
    model.A = sparse([Aeq;Aineq]);
    model.rhs = [beq; bineq];
   
    model.vtype = char(1,3*n+1);
    for i = 1:n+1
        model.vtype(i) = 'I';
    end
    for i = n+2:2*n+1
        model.vtype(i) = 'C';
    end
    for i = 2*n+2:4*n+1
        model.vtype(i) = 'C';
    end
    
    [aux1m, ~] = size(Aeq);
    [aux2m, ~] = size(Aineq);
    model.sense = char(1,aux1m+aux2m);
    for i = 1:aux1m
        model.sense(i) = '=';
    end
    for i = aux1m+1:aux1m+aux2m
        model.sense(i) = '<';
    end
    
    model.modelsense = 'min'; 
    
    result.status = 'INFEASIBLE';
    it2 = 1;
    while strcmp(result.status,'INFEASIBLE')
       result = gurobi(model);
       
       if it2 > nnb
           break;
       end
       
       d = D(problem.intcon,it2);
       Aeq = [problem.Aeq,            zeros(meq,n),           zeros(meq,1), zeros(meq,n),              zeros(meq,n);
            zeros(sd,n),            eyeN(problem.intcon,:), d,             zeros(sd,n),             zeros(sd,n);
            eyeN(problem.intcon,:), -eyeN(problem.intcon,:), zeros(sd,1), eyeN(problem.intcon,:),-eyeN(problem.intcon,:)];
       model.A = sparse([Aeq;Aineq]);
       it2 = it2+1;
    end
       

    if strcmp(result.status,'OPTIMAL')
        f = problem.f'*result.x(1:n);
        fprintf('alg iterations %f\n', result.itercount);
        fprintf('alg function %f\n', result.objval); 
    else
        f = inf;
    end
    
     fprintf('alg exit %s\n', result.status);
     
end