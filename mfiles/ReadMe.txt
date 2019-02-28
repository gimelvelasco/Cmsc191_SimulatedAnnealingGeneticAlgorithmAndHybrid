To Run the Genetic Algorithm, open Velasco_GA.m
To Run the Simulated Annealing, open Velasco_SA.m
To Run the Hybrid Simulated Annealing, open Velasco_HSA.m

On every GA, SA and HSA code found in the mfiles folder, find the block "INPUT ARGUMENTS", the following arguments are found there:

I. Velasco_GA.m:
	= CostF - Objective Function to minimize
		- 1 - DE JONGS
		- 2 - AXIS PARALLEL HYPER-ELLIPSOID
		- 3 - ROTATED HYPER-ELLIPSOID
		- 4 - RASTRIGINS
		- ow - ACKLEYS
	= nVar - number alleles of chromosome
	= VarSize - stores the fittest chromosome per iteration (no need to change)
	= VarMin - lower bound of variable value
	= VarMax - upper bound of variable value
	= MaxGen - maximum number of generations (this is met if tol is not yet met)
	= nPop - chromosome population per generation
	= nSel - Type of Selection process to use
		- 1 - Tournament Selection
		- ow - Roulette Wheel Selection
	= mu - mutation rate
	= sigma - crossover rate
	= cxver - Type of Crossover process to use
		- 1 - Single point Crossover
		- ow - Three point Crossover |
	= tol - fitness tolerance (minimization problem)

II. Velasco_SA.m:
	= CostF - Objective Function to minimize
		- 1 - DE JONGS
		- 2 - AXIS PARALLEL HYPER-ELLIPSOID
		- 3 - ROTATED HYPER-ELLIPSOID
		- 4 - RASTRIGINS
		- ow - ACKLEYS
	= nVar - number alleles of chromosome
	= VarSize - stores the fittest chromosome per iteration (no need to change)
	= VarMin - lower bound of variable value
	= VarMax - upper bound of variable value
	= MaxIt - maximum number of iterations (so as to avoid too much iterations)
	= T0 - initial temperature
	= Tf - final temperature
	= alpha - temperature reduction rate
	= nMove - number of neighbors the solution considers
	= mu - mutation rate

III. Velasco_HSA.m:
	= CostF - Objective Function to minimize
		- 1 - DE JONGS
		- 2 - AXIS PARALLEL HYPER-ELLIPSOID
		- 3 - ROTATED HYPER-ELLIPSOID
		- 4 - RASTRIGINS
		- ow - ACKLEYS
	= nVar - number alleles of chromosome
	= VarSize - stores the fittest chromosome per iteration (no need to change)
	= VarMin - lower bound of variable value
	= VarMax - upper bound of variable value
	= MaxIt - maximum number of iterations (so as to avoid too much iterations)
	= T0 - initial temperature
	= Tf - final temperature
	= alpha - temperature reduction rate
	= nPop - number of chromosomes per iteration
	= nMove - number of neighbors the solution considers
	= mu - mutation rate
	= sigma - crossover rate

The for loop found at the top can be uncommented if the user desires to.
The values found in the arguments are the recommended input arguments for the algorithms.
The raw results of the tests can be found in the latex folder.