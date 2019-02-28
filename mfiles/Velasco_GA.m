%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VELASCO, Gimel David F.
%2012-58922
%Cmsc 191
%Genetic Algorithm
%Final Exam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for runs=1:3
clear;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT ARGUMENTS%%%%%%%Sir Joel, dito po%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CostF = 1; % | 1 - DE JONGS | 2 - AXIS PARALLEL HYPER-ELLIPSOID | 3 - ROTATED HYPER-ELLIPSOID | 4 - RASTRIGINS | ow - ACKLEYS |
nVar = 5; %number of dimensions
VarSize = zeros(nVar); %stores the fittest chromosome per iteration
VarMin = -5.12; %upper bound of variable value
VarMax = 5.12; %lower bound of variable value
MaxGen = 100000; %maximum number of generations
nPop = 5000; %chromosome population
nSel = 1; % | 1 - Tournament Selection | ow - Roulette Wheel Selection |
mu = 0.1; %mutation rate
sigma = 0.25; %crossover rate
cxver= 2; % | 1 - Single point Crossover | ow - Three point Crossover |
tol = 0.0000025; %fitness tolerance (minimization problem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 1: Initialization
test_func = CostF;
ub = VarMax;
lb = VarMin;
ulb = ub; %temporary bengat
tpl = nVar;
popsize = nPop;   %Population Size
maxgens = MaxGen;   %Maximum Generations
pc = sigma;      %crossover rate
pm = mu;       %mutation rate
%tol = tol;  %tolerance
iteration_array = zeros(1);
fittest_array = zeros(1);
solution_array = VarSize;
Chromosome = zeros(popsize,tpl);
for i=1:popsize
    Chromosome(i,:) = 2*ulb*(rand(1,tpl)-0.5);  %initializing first generation
end

%%%%%%%%%%%%%%%%%%%%%%%%%GENETIC ALGORITHM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for generation=1:maxgens   %LOOPS STEPS 2 to 5
%STEP 2: Selection
if nSel == 1
    %Selection by Tournament
    for i=1:popsize/4   %Tournament process
        red = abs(floor(popsize*rand(1))) + 1;  %two random chromosomes to compete
        blu = abs(floor(popsize*rand(1))) + 1;
        if OBJFUNC(Chromosome(red,:),tpl,test_func) > OBJFUNC(Chromosome(blu,:),tpl,test_func)    %competition
            Chromosome(red,:) = Chromosome(blu,:);      %Blue Wins the tournament and overwrites Red
        else
            Chromosome(blu,:) = Chromosome(red,:);      %Red Wins the tournament and overwrites Blue
        end
    end
else
    %Selection by Roulette-Wheel
    F_obj = zeros(1);
    for i=1:popsize
       F_obj(i) = OBJFUNC(Chromosome(i,:),tpl,test_func);
    end
    Total = 0;
    Fitness = zeros(1);
    P = zeros(1);
    C = zeros(1);
    R = zeros(1);
    NewChromosome = zeros(tpl);
    for i=1:popsize
       Fitness(i) = 1/(1+F_obj(i));
       Total = Total + Fitness(i);
    end
    for i=1:popsize
        P(i) = Fitness(i)/Total;
    end
    ctr = 0;
    for i=1:popsize
        ctr = ctr + P(i);
        C(i) = ctr;
    end
    for i=1:popsize
        R(i) = rand(1);
    end
    for i=1:popsize
        NewChromosome(i,:) = zeros(1,tpl);
    for j=1:popsize-1
        if R(i) > C(j) && R(i) <= C(j+1)
            NewChromosome(i,:) = Chromosome(j,:);
        end
    end
    if NewChromosome(i,:) == zeros(1,tpl)
        NewChromosome(i,:) = Chromosome(1,:);
    end
    end
    for i=1:popsize
        Chromosome(i,:) = NewChromosome(i,:);
    end
end
%%
R = rand(1,popsize);
k = zeros(1);
PChromosome = zeros(tpl);
cp = zeros(1);
cp3 = zeros(3);
ctr = 0;    %holds number of parents
for i=1:popsize
    if R(i) < pc
        %select parent
        ctr = ctr + 1;
        k(ctr) = i; %will save the positions of the parent chromosomes
        PChromosome(ctr,:) = Chromosome(i,:);
    end
end
if ctr == 0 %if no parents were selected for the next generation
    continue;
end
%%
%STEP 3: Cross-Over
if cxver == 1 %single point cxver
    for i=1:ctr
        cp(i) = floor(abs((tpl-1)*rand(1)))+1;   %crossover points
    end
    for i=1:ctr-1
       Chromosome(k(i),:) = CROSSOVER(PChromosome(i,:),PChromosome(i+1,:),cp(i),tpl);   %crossover ci and ci+1
    end
    Chromosome(k(ctr),:) = CROSSOVER(PChromosome(ctr,:),PChromosome(1,:),cp(ctr),tpl);    %crossover ck and c1
else %three point cxver
    for i=1:ctr
        cp3(i,:) = floor(abs((tpl-1)*rand(1,3)))+1;   %crossover points
    end
    for i=1:ctr-1
       Chromosome(k(i),:) = CROSSOVER3(PChromosome(i,:),PChromosome(i+1,:),cp3(i,:),tpl);   %crossover ci and ci+1
    end
    Chromosome(k(ctr),:) = CROSSOVER3(PChromosome(ctr,:),PChromosome(1,:),cp3(ctr,:),tpl);    %crossover ck and c1
end
%%
%STEP 4: Mutation
%Per Chromosome mutation
mu = round(pm*popsize); %#ofchromosomestomutate = mutationrate*populationsize
for i=1:mu
    cngn = abs(floor(popsize*rand(1))) + 1; %random popsize number
    q = OBJFUNC(Chromosome(cngn,:),tpl,test_func);
    if q < 1
        Chromosome(cngn,:) = 2*ulb*q*(rand(1,tpl)-0.5);   %mutation
    else
        Chromosome(cngn,:) = 2*ulb*(rand(1,tpl)-0.5);
    end
end
%%
%STEP 5: Post-Evaluation
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
%fprintf('Fittest: %.16f\t\tRuntime: %.2f seconds\n',fittest,toc);
%disp(Chromosome(fi,:))
fittest_array(generation) = fittest;
iteration_array(generation) = generation;
solution_array(generation,:) = Chromosome(fi,:);
if fittest < tol %&& abs(DECODE(Chromosome(fi,:),tpl) - round(DECODE(Chromosome(fi,:),tpl))) < 0.00005
    break;
end
%Step 6: Repeat Generation Iteration
end
%%
%STEP 7: Solution (Best Chromosomes)
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
fprintf('====================GENETIC ALGORITHM=============================\n');
if nSel == 1
fprintf('Tournament Selection. ');
else
fprintf('Roulette-Wheel Selection. ');
end
if cxver == 1
fprintf('Single Point Crossover.\n');
else
fprintf('Triple Point Crossover.\n');
end
fprintf('===================Generations: %d\tPopulation: %d===========\n',maxgens,popsize);
%fprintf('Final Generation\n');
%Chromosome = SORT(Chromosome,popsize,tpl,test_func);
%disp(Chromosome)
fprintf('Generation Number: %d\nPopulation Size: %d\nCrossover Rate: %.2f\nMutation Rate: %.2f\nDimensions: %d\n',maxgens,popsize,pc,pm,tpl);
%fprintf('\nThe Fittest Chromosome is\n');
%disp(Chromosome(fi,:)')
%fprintf('Using Tournament Selection\n');
fprintf('Fitness Function Value of %.16f\nTotal Runtime of %f seconds\n',fittest,toc);
fprintf('The Root for Test Function %d is:\n',test_func);
disp(Chromosome(fi,:))
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
%%
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%