function [i, j, sino_index, accepted_lors] = sinoIndex(options, ring_pos1, ring_pos2, ring_number1, ring_number2)
%SINOINDEX This function obtains the sinogram indices from the input
%detector numbers
%

if options.det_w_pseudo > options.det_per_ring
    ring_pos1 = idivide(ring_pos1, uint16(options.cryst_per_block)) + ring_pos1;
    ring_pos2 = idivide(ring_pos2, uint16(options.cryst_per_block)) + ring_pos2;
end
xa = int16(max([ring_pos1 ring_pos2],[],2));
ya = int16(min([ring_pos1 ring_pos2],[],2));
j = idivide(mod((xa + ya + options.det_w_pseudo/2 + 1), options.det_w_pseudo), 2);
b = j + options.det_w_pseudo/2;
i = abs(xa - ya - options.det_w_pseudo/2);
ind = ya < j | b < xa;
i(ind) = -i(ind);
swap = (j*2) < -i | i <= ((j - options.det_w_pseudo / 2) * 2);
if mod(options.Ndist,2) == 0
    accepted_lors = (i <= (options.Ndist/2 + min(0,options.ndist_side)) & i >= (-options.Ndist/2 + max(0,options.ndist_side)));
else
    accepted_lors = (i <= options.Ndist/2 & i >= (-options.Ndist/2));
end
% accepted_lors = (i < 0 & options.Ndist/2 > (-i - 1)) | (i >= 0 & options.Ndist/2 > i);
j = idivide(j,options.det_w_pseudo/2/options.Nang);
if sum(options.pseudot) > 0
    gapSize = floor(options.rings / (sum(options.pseudot) + 1));
    ring_number1 = idivide(ring_number1, gapSize) + ring_number1;
    ring_number2 = idivide(ring_number2, gapSize) + ring_number2;
end
accepted_lors = accepted_lors & (abs(int16(ring_number1) - int16(ring_number2)) <= options.ring_difference);
ring_number3 = ring_number1;
ring_number1(swap) = ring_number2(swap);
ring_number2(swap) = ring_number3(swap);
clear ring_number3
i = i(accepted_lors);
j = j(accepted_lors);
% swap = swap(accepted_lors);
if min(i) <= 0
    i = i + abs(min(i)) + 1;
end
j = j + 1;
ring_number1 = ring_number1(accepted_lors);
ring_number2 = ring_number2(accepted_lors);
ring_pos1 = ring_pos1(accepted_lors);
ring_pos2 = ring_pos2(accepted_lors);
swappi = ring_pos2 > ring_pos1;
ring_number3 = ring_number1;
ring_number1(swappi) = ring_number2(swappi);
ring_number2(swappi) = ring_number3(swappi);
if options.span == 1
    sino_index = ring_number2 * options.rings + ring_number1;
else
    ind2 = abs(int16(ring_number1) - int16(ring_number2)) <= floor(options.span / 2);
    sino_index = zeros(size(ring_number1,1),1,'int16');
    iD = int16(ring_number1(ind2)) + int16(ring_number2(ind2));
    sino_index(ind2) = iD;
    ind3 = ~ind2;
    erotus = int16(ring_number1(ind3)) - int16(ring_number2(ind3));
    summa = int16(ring_number1(ind3)) + int16(ring_number2(ind3));
    index = idivide(abs(erotus) + floor(options.span / 2), int16(options.span));
    seg = cumsum(options.segment_table);
    ind4 = erotus < 0;
    index1 = (summa(ind4) - (floor(options.span/2) .* (index(ind4) * 2 - 1) + index(ind4))) + int16(seg((index(ind4)) * 2 - 1)');
    ind4 = erotus > 0;
    index2 = (summa(ind4) - (floor(options.span/2) .* (index(ind4) * 2 - 1) + index(ind4))) + int16(seg((index(ind4)) * 2)');
    index_ = zeros(length(ind4),1,'int16');
    index_(erotus < 0) = index1;
    index_(erotus > 0) = index2;
    sino_index(ind3) = index_;
end
sino_index = sino_index + 1;
end

