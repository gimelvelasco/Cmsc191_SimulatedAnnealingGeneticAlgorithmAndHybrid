%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Simulated Annealing
%Final Exam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for runs=1:3
clear;
tic;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT ARGUMENTS%%%%%%%Sir Joel, dito po%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CostF = 1; % | 1 - DE JONGS | 2 - AXIS PARALLEL HYPER-ELLIPSOID | 3 - ROTATED HYPER-ELLIPSOID | 4 - RASTRIGINS | ow - ACKLEYS |
nVar = 5; %number of dimensions
VarSize = zeros(nVar);
VarMin = -5.12; %upper bound of variable value
VarMax = 5.12; %lower bound of variable value
MaxIt = 100000;
T0 = 100;
Tf = 0.000000000000001;
alpha = 0.7;
%nPop = 10000;
nMove = 10000;
mu = 0.1;
%sigma = 0.25;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
test_func = CostF;  %sets the number of w/c test function to be solved
ub = VarMax;
lb = VarMin;
ulb = ub;        %upper and lower bound
tpl = nVar;      %dimensions
x_sol = 2*ulb*(rand(1,tpl)-0.5);               %initial guess of the solution x
cooling_ratio = alpha;                    %sets the cooling ratio to 0.8  i.e. 0.7 < 0.8 < 0.9
num_neigh = nMove;                      %initializes the size of the random neighbors
cooling_sched = zeros(1);               %pre-allocation for speed
cooling_sched(1) = T0;                 %initializes the cooling schedule T0
iteration_array = zeros(1);
fittest_array = zeros(1);
solution_array = VarSize;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%SIMULATED ANNEALING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sched = 1;                                  %index
while cooling_sched(sched) > Tf     %iteration will stop if the cooling temperature reached less than 0.00000001
    T = cooling_sched(sched);               %sets the value of the temperature T
    for j=1:num_neigh
        r  = (cooling_ratio)^sched;             %is used so that the randomness of selecting a neighbor becomes narrower
        x_tmp = 2*ulb*r*(rand(1,tpl)-0.5);    %selects a random neighbor for comparison. with a decreasing amount of randomness
        if OBJFUNC(x_tmp,tpl,test_func) < OBJFUNC(x_sol,tpl,test_func)  %if the neighbor is better, change the solution
            x_sol = x_tmp;
        elseif OBJFUNC(x_tmp,tpl,test_func) > OBJFUNC(x_sol,tpl,test_func)  %if not, change the solution if it is lucky
            delta = OBJFUNC(x_tmp,tpl,test_func) - OBJFUNC(x_sol,tpl,test_func);
            p = P(delta,T);
            q = rand(1);
            if q <= p
                x_sol = x_tmp; 
            end
        end
    end

	fittest_array(sched) = OBJFUNC(x_sol,tpl,test_func);
	iteration_array(sched) = sched;
	solution_array(sched,:) = x_sol;

    cooling_sched(sched+1) = T*(cooling_ratio)^sched;
    sched = sched+1;
    if sched > MaxIt
    	break;
    end
end
%%
%SOLUTION
if test_func == 1
fprintf('====================DE JONGS FUNCTION=============================\n');
elseif test_func == 2
fprintf('==============AXIS PARALLEL HYPER-ELLIPSOID FUNCTION==============\n');
elseif test_func == 3
fprintf('===============ROTATED HYPER-ELLIPSOID FUNCTION====================\n');
elseif test_func == 4
fprintf('====================RASTRIGINS FUNCTION===========================\n');
else
fprintf('=====================ACKLEYS FUNCTION=============================\n');
end
fprintf('==================SIMULATED ANNEALING=============================\n');
fprintf('With the Objective Function Value of %.16f\nTotal Runtime of %f seconds\nAnd Final Cooling Temperature of %.16f\n',OBJFUNC(x_sol,tpl,test_func),toc,cooling_sched(sched));
fprintf('The Root for Test Function %d is\n',test_func);
disp(x_sol)
fprintf('==================================================================\n');
%%
figure
subplot(2,1,1);
plot(iteration_array,fittest_array);
legend('Cost Function Value');
xlabel('Generation');
ylabel('Fitness of fittest Chromosome');

solution_array = transpose(solution_array);
subplot(2,1,2);
plot(iteration_array,solution_array);
xlabel('Generation');
ylabel('Solution of fittest Chromosome');
%end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%