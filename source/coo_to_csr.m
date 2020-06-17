function csr_row = coo_to_csr(coo_row)
% This function converts sparse COO-format row indices into sparse
% CSR-format row indices. Output is in unsigned integer format with the
% same byte count as input data, (e.g. single in unsigned 32-bit).

[~,IA] = ismember(coo_row,coo_row);
IA = unique(IA);
if isa(coo_row,'int8') || isa(coo_row,'uint8')
    csr_row = [uint8(IA-1);uint8(length(coo_row))];
elseif isa(coo_row,'int16') || isa(coo_row,'uint16')
    csr_row = [uint16(IA-1);uint16(length(coo_row))];
elseif isa(coo_row,'uint32') || isa(coo_row,'int32') || isa(coo_row,'single')
    csr_row = [uint32(IA-1);uint32(length(coo_row))];
else
    csr_row = [uint64(IA-1);uint64(length(coo_row))];
end
end