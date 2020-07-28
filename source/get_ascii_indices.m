function ascii_ind = get_ascii_indices(coincidence_mask)
%GET_ASCII_INDICES Get the necessary column indices for the ASCII files
%   Converts the GATE ASCII coincidence mask (setCoincidenceMask) to
%   correspond to the column indices of the data.

if (length(coincidence_mask) < 36 || length(coincidence_mask) > 36) && ~isempty(coincidence_mask)
    error('Length of the coincidence mask is too little or too much. Needs to be exactly 36 or an empty array')
end
if isempty(coincidence_mask)
    coincidence_mask = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
end
ascii_ind.run_index1 = 0;
ascii_ind.run_index2 = 0;
ascii_ind.event_index1 = 0;
ascii_ind.event_index2 = 0;
ascii_ind.source_index1 = 0;
ascii_ind.source_index2 = 0;
ascii_ind.time_index = 0;
ascii_ind.time_index2 = 0;
ascii_ind.energy_index1 = 0;
ascii_ind.world_index1 = 0;
ascii_ind.base_ind1 = 0;
ascii_ind.rsector_ind1 = 0;
ascii_ind.module_ind1 = 0;
ascii_ind.submodule_ind1 = 0;
ascii_ind.crs_ind1 = 0;
ascii_ind.layer_ind1 = 0;
ascii_ind.scatter_index_cp1 = 0;
ascii_ind.scatter_index_cd1 = 0;
ascii_ind.scatter_index_rp1 = 0;
ascii_ind.scatter_index_rd1 = 0;
ascii_ind.scanner_axial_index1 = 0;
ascii_ind.scanner_angular_index1 = 0;
ascii_ind.energy_index2 = 0;
ascii_ind.world_index2 = 0;
ascii_ind.base_ind2 = 0;
ascii_ind.rsector_ind2 = 0;
ascii_ind.module_ind2 = 0;
ascii_ind.submodule_ind2 = 0;
ascii_ind.crs_ind2 = 0;
ascii_ind.layer_ind2 = 0;
ascii_ind.scatter_index_cp2 = 0;
ascii_ind.scatter_index_cd2 = 0;
ascii_ind.scatter_index_rp2 = 0;
ascii_ind.scatter_index_rd2 = 0;
ascii_ind.scanner_axial_index2 = 0;
ascii_ind.scanner_angular_index2 = 0;
summa = 0;
pituus = length(coincidence_mask);
for kk = 1 : pituus
   if kk == 1 && coincidence_mask(kk) == 1
       ascii_ind.run_index1 = 1;
       summa = summa + 1;
   elseif kk == 2 && coincidence_mask(kk) == 1
       ascii_ind.event_index1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 3 && coincidence_mask(kk) == 1
       ascii_ind.sourceid_index1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 4 && coincidence_mask(kk) == 1
       ascii_ind.source_index1 = 1 + summa;
       summa = summa + 3;
   elseif kk == 5 && coincidence_mask(kk) == 1 && coincidence_mask(kk-1) == 0
       error('Only one source index is selected! All three need to be selected')
   elseif kk == 6 && coincidence_mask(kk) == 1 && (coincidence_mask(kk-1) == 0 || coincidence_mask(kk-2) == 0)
       error('Only one source index is selected! All three need to be selected')
   elseif kk == 7 && coincidence_mask(kk) == 1
       ascii_ind.time_index = 1 + summa;
       summa = summa + 1;
   elseif kk == 8 && coincidence_mask(kk) == 1
       ascii_ind.energy_index1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 9 && coincidence_mask(kk) == 1
       ascii_ind.world_index1 = 1 + summa;
       summa = summa + 3;
   elseif kk == 12 && coincidence_mask(kk) == 1
       ascii_ind.base_ind1 = 1 + summa;
       summa = summa + 1;
       ascii_ind.rsector_ind1 = 1 + summa;
       summa = summa + 1;
       ascii_ind.module_ind1 = 1 + summa;
       summa = summa + 1;
       ascii_ind.submodule_ind1 = 1 + summa;
       summa = summa + 1;
       ascii_ind.crs_ind1 = 1 + summa;
       summa = summa + 1;
       ascii_ind.layer_ind1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 13 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_cp1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 14 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_cd1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 15 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_rp1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 16 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_rd1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 17 && coincidence_mask(kk) == 1
       ascii_ind.scanner_axial_index1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 18 && coincidence_mask(kk) == 1
       ascii_ind.scanner_angular_index1 = 1 + summa;
       summa = summa + 1;
   elseif kk == 19 && coincidence_mask(kk) == 1
       ascii_ind.run_index2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 20 && coincidence_mask(kk) == 1
       ascii_ind.event_index2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 21 && coincidence_mask(kk) == 1
       ascii_ind.sourceid_index2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 22 && coincidence_mask(kk) == 1
       ascii_ind.source_index2 = 1 + summa;
       summa = summa + 3;
   elseif kk == 23 && coincidence_mask(kk) == 1 && coincidence_mask(kk-1) == 0
       error('Only one source index is selected! All three need to be selected')
   elseif kk == 24 && coincidence_mask(kk) == 1 && (coincidence_mask(kk-1) == 0 || coincidence_mask(kk-2) == 0)
       error('Only one source index is selected! All three need to be selected')
   elseif kk == 25 && coincidence_mask(kk) == 1
       ascii_ind.time_index2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 26 && coincidence_mask(kk) == 1
       ascii_ind.energy_index2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 27 && coincidence_mask(kk) == 1
       ascii_ind.world_index2 = 1 + summa;
       summa = summa + 3;
   elseif kk == 30 && coincidence_mask(kk) == 1
       ascii_ind.base_ind2 = 1 + summa;
       summa = summa + 1;
       ascii_ind.rsector_ind2 = 1 + summa;
       summa = summa + 1;
       ascii_ind.module_ind2 = 1 + summa;
       summa = summa + 1;
       ascii_ind.submodule_ind2 = 1 + summa;
       summa = summa + 1;
       ascii_ind.crs_ind2 = 1 + summa;
       summa = summa + 1;
       ascii_ind.layer_ind2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 31 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_cp2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 32 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_cd2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 33 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_rp2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 34 && coincidence_mask(kk) == 1
       ascii_ind.scatter_index_rd2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 35 && coincidence_mask(kk) == 1
       ascii_ind.scanner_axial_index2 = 1 + summa;
       summa = summa + 1;
   elseif kk == 36 && coincidence_mask(kk) == 1
       ascii_ind.scanner_angular_index2 = 1 + summa;
       summa = summa + 1;
   end
end
end

