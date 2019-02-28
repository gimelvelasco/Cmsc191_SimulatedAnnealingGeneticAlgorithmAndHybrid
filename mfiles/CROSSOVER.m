function Chromosome_ret = CROSSOVER(PChromosome1,PChromosome2,cp_rec,tpl)
    Chromosome_ret = PChromosome1;
    for i=cp_rec+1:tpl    %crossovers after the crossover point
        Chromosome_ret(i) = PChromosome2(i);
    end
end