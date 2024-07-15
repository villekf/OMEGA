function output = ECOSEM(im, D, OSEM_apu, COSEM_apu, epps)
alpha_eco = 1;
output = alpha_eco .* OSEM_apu + (1 - alpha_eco) .* COSEM_apu;
eco_s1 = sum(D .* (-COSEM_apu .* log(im + epps) + im));
eco_s2 = sum(D .* (-COSEM_apu .* log(output + epps) + output));
while (alpha_eco > 0.0096 && eco_s1 < eco_s2)
    alpha_eco = alpha_eco .* 0.9;
    output = alpha_eco .* OSEM_apu + (1 - alpha_eco) .* COSEM_apu;
    eco_s2 = sum(D .* (-COSEM_apu .* log(output + epps) + output));
end
if (alpha_eco <= 0.0096)
    output = COSEM_apu;
end
end