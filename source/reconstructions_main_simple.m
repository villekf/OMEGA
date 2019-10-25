function f_osem = reconstructions_main_simple(options)
%reconstructions_main_simple Reconsruction file for simple reconstruction
%   Extract OSEM data
pz = reconstructions_main(options);
f_osem = pz{2};
end