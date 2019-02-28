function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec,tpl)
	%bubblesort cp_rec
	for i=1:2
		for j=1:2
			if cp_rec(j) > cp_rec(j+1)
				tmp = cp_rec(j);
				cp_rec(j) = cp_rec(j+1);
				cp_rec(j+1) = tmp;
			end
		end
	end

	%crossover process
    Chromosome_ret = PChromosome1;
    for i=cp_rec(1)+1:cp_rec(2)+1    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
    for i=cp_rec(2)+1:cp_rec(3)+1    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome1(i);
    end
    for i=cp_rec(3)+1:tpl
        Chromosome_ret(i) = PChromosome2(i);
    end
end