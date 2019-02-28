function bp_ret = BOLTZMAN(Ei,Ej,T)
    bp_ret = 1/(1 + exp((Ei-Ej)/T));
end