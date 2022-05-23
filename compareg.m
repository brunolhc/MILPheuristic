function[] = compareg()
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab\x64_win64')
addpath('C:\gurobi910\win64\matlab')

it = 150;

nome = strcat('results.csv');
arquivo = fopen(nome, 'w+');

fprintf(arquivo, 'Problem m n f_alg t_alg it_alg f_g t_g it_g');
fprintf(arquivo, '\n');

for i = 1:it
    
     prob = {'App2-1.mps';
'sp98ir.mps';
'haprp.mps';
'enlight_hard.mps';
'App2-2.mps';
'mod010.mps';
'RococoC10-001000.mps';
'f2gap201600.mps';
'f2gap401600.mps';
'f2gap801600.mps';
'Eil33-2.mps';
'gt2.mps';
'p0201.mps';
'supportcase14.mps';
'supportcase16.mps';
'Chromaticindex32-8.mps';
'Neos-1516309.mps';
'Neos-1599274.mps';
'Neos-831188.mps';
'neos-3592146-hawea.mps';
'Neos-555424.mps';
'neos-3381206-awhea.mps';
'Neos-555343.mps';
'Neos-555001.mps';
'neos-2624317-amur.mps';
'neos18.mps';
'manna81.mps';
'mzzv11.mps';
'enlight8.mps';
'f2gap40400.mps';
'pw-myciel4.mps';
'acc-tight5.mps';
'opm2-z6-s1.mps';
'acc-tight2.mps';
'acc-tight4.mps';
'opm2-z7-s8.mps';
'p2756.mps';
'harp2.mps';
'disktom.mps'
};

    fprintf('\n Problema %s \n', char(prob(i)));
    %name = strcat('problem',num2str(i),'.mat');
    problem = mpsread(char(prob(i)));
    
    [m,n] = size([problem.Aineq;problem.Aeq]);
    [meq, ~] = size(problem.Aeq);
    [mineq, ~] = size(problem.Aineq);
    
    clear model;
    model.A = sparse([problem.Aeq;problem.Aineq]);
    model.obj = problem.f;
    model.rhs = [problem.beq;problem.bineq];
    model.vtype = char(1,n);
    for j = 1:n
        model.vtype(j) = 'I';
    end

    model.sense = char(1,meq+mineq);
    for j = 1:meq
        model.sense(j) = '=';
    end
    for j = meq+1:meq+mineq
        model.sense(j) = '<';
    end

    model.ub = problem.ub;
    model.lb = problem.lb;
    model.modelsense = 'min';
    
    tic
    result = gurobi(model);
    t = toc;
    
    logcplexnome = strcat(char(prob(i)),'-gurobi.txt');
    loggurobi = fopen(logcplexnome, 'w+');
    fprintf(loggurobi, 'gurobi iterations %f\n', result.itercount);
    fprintf(loggurobi, 'gurobi function %f\n', result.objval);
    fprintf(loggurobi, 'gurobi exit %s\n \n', result.status);
    
    fprintf(loggurobi, 'size Aeq: %d\n \n', size(problem.Aeq));
    
    if strcmp(result.status,'OPTIMAL')
        pos = abs(result.x - round(result.x)) > 1e-4;
        nint = result.x(pos);
        aux = 1:n;
        npos = aux(pos);
        fprintf(loggurobi, 'x = [');
        for j = 1:n
            fprintf(loggurobi, '%f ', result.x(j));
        end
        fprintf(loggurobi, '] \n \n');
        
        fprintf(loggurobi, 'non integer: \n');
        for j = 1:length(npos)
            fprintf(loggurobi, 'pos: %d, val:%f \n', npos(j), nint(j));
        end
    end
    
    fprintf('gurobi iterations %f\n', result.itercount);
    fprintf('gurobi function %f \n', result.objval);
    fprintf('gurobi exit %s\n \n', result.status);
    
    tic
    [resultalg, f] = alg2g(problem);
    t_alg = toc;
    
    logcplexnome = strcat(char(prob(i)),'-alg.txt');
    logalg = fopen(logcplexnome, 'w+');
    fprintf(logalg, 'alg iterations %f\n', resultalg.itercount);
    fprintf(logalg, 'alg function %f\n', f);
    fprintf(logalg, 'alg exit %s\n \n', resultalg.status);
    
    fprintf(logalg, 'size Aeq: %d\n \n', size(problem.Aeq));
    
    if strcmp(resultalg.status,'OPTIMAL')
        x = resultalg.x(1:n);
        pos = abs(x - round(x)) > 1e-4;
        nint = x(pos);
        aux = 1:n;
        npos = aux(pos);
        fprintf(logalg, 'x = [');
        for j = 1:n
            fprintf(logalg, '%f ', x(j));
        end
        fprintf(logalg, '] \n \n');
        
        fprintf(logalg, 'non integer: ');
        z = 1;
        for j = 1:length(npos)
            z= 0;
            fprintf(logalg, 'pos: %d, val:%f \n', npos(j), nint(j));
        end
        if z
            fprintf(logalg, 'FALSE \n');
        end
            
        breq = problem.Aeq*resultalg.x(1:n);
        brineq = problem.Aineq*resultalg.x(1:n);
        
        for k = 1:meq
            if breq(k) < problem.beq(k) - 1e-4
                fprintf(logalg, '\n Aeq violation row %d: ax = %f; beq = %f', k, breq(k), problem.beq(k));
            end
            if breq(k) > problem.beq(k) + 1e-4
                fprintf(logalg, '\n Aeq violation row %d: ax = %f; beq = %f', k, breq(k), problem.beq(k));
            end
        end
        
        for k = 1:mineq
            if brineq(k) > problem.bineq(k) + 1e-4
                fprintf(logalg, '\n Aineq violation row %d: ax = %f; bineq = %f', k, brineq(k), problem.bineq(k));
            end
        end
        
        for k = 1:n
            if x(k) > problem.ub(k) + 1e-4
                fprintf(logalg, '\n ub violation row %d: x = %f; ub = %f', k, x(k), problem.ub(k));
            end
            if x(k) < problem.lb(k) - 1e-4
                fprintf(logalg, '\n lb violation row %d: x = %f; lb = %f', k, x(k), problem.lb(k));
            end
        end
        

        
%         if problem.Aeq*resultalg.x(1:n) >= problem.beq - 1e-4
%             fprintf(logalg, '\n Aeq*x >= beq: TRUE \n');
%         else
%             fprintf(logalg, '\n Aeq*x >= beq: FALSE \n');
%         end
%         if problem.Aeq*resultalg.x(1:n) <= problem.beq + 1e-4
%             fprintf(logalg, '\n Aeq*x <= beq: TRUE \n');
%         else
%             fprintf(logalg, '\n Aeq*x <= beq: FALSE \n');
%         end
%         
%         if problem.Aineq*resultalg.x(1:n) <= problem.bineq + 1e-4
%             fprintf(logalg, '\n Aineq*x <= bineq: TRUE \n');
%         else 
%             fprintf(logalg, '\n Aineq*x <= bineq: FALSE \n');
%         end
%         
%         if resultalg.x(1:n) <= problem.ub
%             fprintf(logalg, '\n x <= ub: TRUE \n');
%         else 
%             fprintf(logalg, '\n x <= ub: FALSE \n');
%         end
%         
%         if resultalg.x(1:n) >= problem.lb
%             fprintf(logalg, '\n x >= lb: TRUE \n');
%         else 
%             fprintf(logalg, '\n x >= lb: FALSE \n');
%         end
    end
    
    if f == 0
        it_cplex = 0;
    else
        it_cplex = result.itercount;
    end
    
    fprintf(arquivo, '%s %d %d %.1f %.1f %d %.1f %.1f %d', char(prob(i)), m, n, f, t_alg, resultalg.itercount, result.objval, t, it_cplex);
    fprintf(arquivo, '\n');
    
    fclose(logalg);
    fclose(loggurobi);
    
end
fclose(arquivo);