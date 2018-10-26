function SinM = form_sinograms(options)
%% FORM MICHELOGRAMS AND THEN SINOGRAMS FROM RAW DETECTOR PAIR DATA
% This code forms the sinograms for the current machine from the raw
% detector pairs. It first creates the Michelogram and later the sinograms
% with the dimensions specified by the user.
% Input the machine and sinogram properties and the raw detector pair data
% Output is the sinograms for each time point

machine_name = options.machine_name;
name = options.name;
rings = options.rings;
Nang = options.Nang;
Ndist = options.Ndist;
pseudot = options.pseudot;
segment_table = options.segment_table;
span = options.span;
NSlices = options.TotSinos;
Nz = options.Nz;
ring_difference = options.ring_difference;
partitions = options.partitions;
tot_time = options.tot_time;

ringsp = rings + length(pseudot);

if exist( 'options.coincidences' , 'var') == 0
    if options.partitions == 1
        if options.use_ASCII
            load([machine_name '_measurements_' name '_static_ASCII.mat'], 'coincidences')
        elseif options.use_LMF
            load([machine_name '_measurements_' name '_static_LMF.mat'], 'coincidences')
        elseif options.use_root
            load([machine_name '_measurements_' name '_static_root.mat'], 'coincidences')
        end
    else
        if options.use_ASCII
            load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_ASCII.mat'], 'coincidences')
        elseif options.use_LMF
            load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_LMF.mat'], 'coincidences')
        elseif options.use_root
            load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_root.mat'], 'coincidences')
        end
    end
end

if options.verbose
    tic
end

data = load([machine_name '_app_coordinates_' num2str(Ndist) 'x' num2str(Nang) '.mat']);

i = data.i;
j = data.j;
joku12 = data.joku12;
clear data


SinM = cell(partitions,1);

for llo = 1 : partitions
    
    P1 = options.coincidences{llo};
    
    Sinog = cell(ringsp,ringsp);
    
    ix=cellfun('isempty',Sinog);
    Sinog(ix)={zeros(Nang,Ndist,'uint16')};
    
    % Create the Michelograms
    tic
    for ii=1:ringsp
        if any(ii==pseudot)
            continue
        end
        for jj=1:ringsp
            if any(jj==pseudot)
                continue
            end
            if issparse(P1{ii,jj})
                CC = uint16(full(P1{ii,jj}));
            else
                CC = P1{ii,jj};
            end
            CC = CC(joku12);
            Sinog{ii,jj} = uint16(accumarray([j i],CC,[Nang Ndist]));
        end
        disp([num2str(100*ii/ringsp) '%'])
    end
    toc
    
    %%
    Sinog = cat(3,Sinog{:});
    Sin = zeros(Nang,Ndist,NSlices,'uint16');
    
    
    
    kkj = [];
    for kk = 0 : floor((ring_difference-ceil(span/2))/span)
        kkj = [kkj; ceil(span/2) + span*kk];
    end
    offset2 = cumsum(segment_table);
    % Create the sinograms
    % First the detectors on the same ring
    Sin(:,:,1:2:Nz) = Sinog(:,:,1:ringsp+1:ringsp^2);
    % Then the detectors on adjecent rings
    for j=1:floor(span/2)
        apu = Sinog(:,:,j*ringsp+1:ringsp+1:ringsp^2);
        apu2 = Sinog(:,:,j+1:ringsp+1:(ringsp-j)*ringsp);
        Sin(:,:,j+1:2:offset2(1)-j) = Sin(:,:,j+1:2:offset2(1)-j) + apu + apu2;
    end
    % Lastly the rest of the detectors with the amount of combined LORs
    % specified with the span value
    for i=1:floor(length(segment_table)/2)
        for j=1:span
            apu = Sinog(:,:,(kkj(i)+j-1)*ringsp+1:ringsp+1:end);
            Sin(:,:,offset2(2*i-1)+j:2:offset2(2*i)-j+1) = Sin(:,:,offset2(2*i-1)+j:2:offset2(2*i)-j+1) + (apu);
            apu2 = Sinog(:,:,kkj(i)+j:ringsp+1:(ringsp-kkj(i)-j+1)*ringsp);
            Sin(:,:,offset2(2*i)+j:2:offset2(2*i+1)-j+1) = Sin(:,:,offset2(2*i)+j:2:offset2(2*i + 1)-j+1) + (apu2);
        end
    end
    SinM{llo} = Sin;
end
%%
if partitions == 1
    save([machine_name '_' name '_sinograms_combined_static_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '.mat'],'SinM')
else
    save([machine_name '_' name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '.mat'], 'SinM')
end
if options.verbose
    disp('Sinograms formed')
    toc
end
end