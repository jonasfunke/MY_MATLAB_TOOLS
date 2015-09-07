function [  ] = export_to_igor_cell( cell_data, wave_names, file_location )
% Writes an txt file that can be imported to igor pro

    if length(cell_data)== length(wave_names)
        % write header/wave names
        file=fopen(file_location, 'w'); %open file to write
        for j=1:length(cell_data)
            fprintf(file, [wave_names{j} '\n']);
            for i=1:length(cell_data{j})                                 %write wavenames at each column header
                fprintf(file, '%e\n',cell_data{j}(i));
            end
            fprintf(file,'\n');
        end
        fclose(file);

        % append data
%         dlmwrite(file_location, data, 'delimiter', '\t','-append')
    else
        disp('Error: columns not equals number of wave_names.')
    end

end
