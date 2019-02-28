%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Hybrid Simulated Annealing
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
nVar = 5;
VarSize = zeros(nVar);
VarMin = -5.12; %upper bound of variable value
VarMax = 5.12; %lower bound of variable value
MaxIt = 100000;
T0 = 100000;
Tf = 0.0000000000000001;
alpha = 0.7;
nPop = 1000;
nMove = 50;
mu = 0.05;
sigma = 0.9;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Initialization Section%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_T = T0;                       %initial temperature of the system
cooling_stop = Tf;      %cooling stops when it reaches this temperature
test_func = CostF;  %sets the number of w/c test function to be solved
popsize = nPop;                         %Population Size
pc = sigma;                               %crossover rate
pm = mu;                               %mutation rate
cooling_ratio = alpha;      %sets the cooling ratio to 0.8  i.e. 0.7 < 0.8 < 0.9
ub = VarMax;
lb = VarMin;
ulb = ub;        %upper and lower bound
tpl = nVar;        %dimensions
num_neigh = nMove;                 %number of neighbors to consider
iteration_array = zeros(1);
fittest_array = zeros(1);
solution_array = VarSize;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%HYBRID SIMULATED ANNEALING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I
cooling_sched = zeros(1);               %pre-allocation for speed
cooling_sched(1) = initial_T;                 %initializes the cooling schedule T0
%II
Chromosome = zeros(popsize,tpl);
for i=1:popsize
    Chromosome(i,:) = 2*ulb*(rand(1,tpl)-0.5);  %initializing first generation
end
%%
sched = 1;                              %index / iteration
while cooling_sched(sched) > cooling_stop     %iteration will stop if the cooling temperature reached less than 0.00000001
    T = cooling_sched(sched);               %sets the value of the temperature T
    %III.a. Do N/2 times
    for j=1:popsize/2
        %III.a.i. Select two parents at random
        red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
        blu = abs(floor(popsize*rand(1))) + 1;  %red and blu hold the index of the parents
        %III.a.ii. Generate two offsprings
        %%Recombination Operator (CROSSOVER)
        pc_trial = rand(1);
        if pc_trial < pc     %if trial made it in the crossover rate
            cp = floor(abs((tpl-1)*rand(1)))+1;   %random crossover point
            Child_Chromosome(1,:) = CROSSOVER(Chromosome(red,:),Chromosome(blu,:),cp,tpl);     %crossover red and blu
            Child_Chromosome(2,:) = CROSSOVER(Chromosome(blu,:),Chromosome(red,:),cp,tpl);     %they will have two children
            %%Neighborhood Operator (MUTATION)
            for k=1:2
                x_sol = Child_Chromosome(k,:);               
                for i=1:num_neigh
                    adrs = abs(floor(popsize*rand(1))) + 1; %gets a random address of a neighbor within the population
                    x_tmp = Chromosome(adrs,:);    %selects a random neighbor for comparison. with a decreasing amount of randomness
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
                Child_Chromosome(k,:) = x_sol;           %will overwrite the child based on the neighborhood operator
            end
            %%III.a.iii. Boltzman Trials
            ARpossibility = rand(1);                    % <0.5 - Single Acceptance/Rejection | >=0.5 - Double Acceptance/Rejection
            if ARpossibility < 0.5 %%Case 1: Double Acceptance/Rejection
                E1 = OBJFUNC(Chromosome(red,:),tpl,test_func) + OBJFUNC(Chromosome(blu,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(1,:),tpl,test_func) + OBJFUNC(Child_Chromosome(2,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(red,:) = Child_Chromosome(1,:);
                    Chromosome(blu,:) = Child_Chromosome(2,:);
                end
            else %%Case 2: Single Acceptance/Rejection
                E1 = OBJFUNC(Chromosome(red,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(1,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(red,:) = Child_Chromosome(1,:);  %offsprings wins the trial
                end

                E1 = OBJFUNC(Chromosome(red,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(2,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(red,:) = Child_Chromosome(2,:);  %offsprings wins the trial
                end

                E1 = OBJFUNC(Chromosome(blu,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(1,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(blu,:) = Child_Chromosome(1,:);  %offsprings wins the trial
                end

                E1 = OBJFUNC(Chromosome(blu,:),tpl,test_func);
                E2 = OBJFUNC(Child_Chromosome(2,:),tpl,test_func);
                bp = BOLTZMAN(E1,E2,T);
                bp_trial = rand(1);
                if bp_trial >= bp
                    %%III.a.iv. Overwrite Parents with the Trial Winner
                    Chromosome(blu,:) = Child_Chromosome(2,:);  %offsprings wins the trial
                end
            end
            %%
        else                    %if the whole trial did not make it inside the crossover rate, it will have a tournament
            if OBJFUNC(Chromosome(red,:),tpl,test_func) > OBJFUNC(Chromosome(blu,:),tpl,test_func)    %competition
                Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
            else
                Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
            end
        end
    end
    %III.b. Periodically Lower T
    %%
    %Post-Iteration-Evaluation
    F_obj = zeros(1);
    for i=1:popsize
       F_obj(i) = OBJFUNC(Chromosome(i,:),tpl,test_func);
    end
    %%
    fittest = F_obj(1);
    for i=1:popsize
        fi = 1;
        if fittest > F_obj(i)
           fittest = F_obj(i);
           fi = i;
        end
    end
    %%
    fittest_array(sched) = fittest;
    iteration_array(sched) = sched;
    solution_array(sched,:) = Chromosome(fi,:);

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
fprintf('==================HYBRID SIMULATED ANNEALING=======================\n');
fprintf('Population Size: %d\tCrossover Rate: %.2f\tMutation Rate: %.2f\tCooling Ratio: %.1f\n',popsize,pc,pm,cooling_ratio);
fprintf('Minimum Energy:\t\t\t\t\t\t\t\t\t%.16f\nTotal Runtime:\t\t\t\t\t\t\t\t\t%f seconds\nFinal Temperature:\t\t\t\t\t\t\t\t%.16f\n',fittest,toc,cooling_sched(sched));
fprintf('Global Minimum is at:\t\t\t\t\t\t');
disp(Chromosome(fi,:))
fprintf('===================================================================\n');
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
%%
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%