function options = set_GATE_variables(options)
%set_GATE_variables Set the GATE variables for non-GATE data
%   If some of the fieldnames in the 'options' struct are missing, errors will
%   be thrown. This function will fill those fieldnames in order to prevent
%   these errors.
if ~isfield(options, 'use_ASCII')
    options.use_ASCII = false;
end
if ~isfield(options, 'use_binary')
    options.use_binary = false;
end
% options.rsector_ind1 = 0;
% options.rsector_ind2 = 0;
% options.module_ind1 = 0;
% options.module_ind2 = 0;
% options.crs_ind1 = 0;
% options.crs_ind2 = 0;
% options.time_index = 0;
% options.event_index1 = 0;
% options.event_index2 = 0;
% options.source_index1 = 0;
% options.source_index2 = 0;
% options.scatter_index1 = 0;
% options.scatter_index2 = 0;
if ~isfield(options, 'use_LMF')
    options.use_LMF = false;
end
options.header_bytes = 0;
options.data_bytes = 0;
options.R_bits = 0;
options.M_bits = 0;
options.S_bits = 0;
options.C_bits = 0;
options.L_bits = 0;
options.coincidence_window = 10e-9;
options.clock_time_step = 1e-12;
if ~isfield(options, 'use_root')
    options.use_root = false;
end
if ~isfield(options, 'fpath')
    options.fpath = '';
end
end

