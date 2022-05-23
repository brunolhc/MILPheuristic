function[problem] = problemgen(m, n)

problem.Aineq = zeros(m,n);
problem.f = zeros(n,1);
problem.bineq = zeros(m,1);

for i = 1:m
    problem.bineq(i) = randi(30*m*n);
    for j = 1:n
        problem.Aineq(i,j) = randi(m*n);
        problem.f(j) = -randi(n);
    end
end

problem.lb = zeros(n,1);
problem.ub = inf(1,n);
problem.Aeq = zeros(0,n);
problem.beq = zeros(0,1);
problem.intcon = 1:n;



        