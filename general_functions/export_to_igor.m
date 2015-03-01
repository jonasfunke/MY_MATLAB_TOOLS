function [  ] = export_to_igor( data, wave_names, file_location )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if size(data,2)== length(wave_names)
        % write header/wave names
        file=fopen(file_location, 'w'); %open file to write
        for i=1:size(data,2)                                 %write wavenames at each column header
            fprintf(file, [wave_names{i} '\t']);
        end
        fprintf(file,'\n');
        fclose(file);

        % append data
        dlmwrite(file_location, data, 'delimiter', '\t','-append')
    else
        disp('Error: columns not equals number of wave_names.')
    end

end

