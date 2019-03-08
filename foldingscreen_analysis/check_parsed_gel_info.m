function [parsed] = check_parsed_gel_info(parsed_data)
    %% check if all important objects are there
    parsed = false;
    fields_required = {'user', 'project', 'design_name', 'date', 'scaffold_type', ...
        'scaffold_concentration', 'staple_concentration', 'comment', ...
        'lanes', 'gelsize', 'agarose_concentration', 'staining', 'mg_concentration', ...
        'voltage', 'running_time', 'cooling'};

    missing_fields = {};

    for i=1:length(fields_required)
        if ~isfield(parsed_data, fields_required{i})
            disp(['Warning: ' fields_required{i} ' not found. Check the gel_info.txt file.'])
            missing_fields = [missing_fields fields_required{i}];

        end
    end

    if isempty(missing_fields)
        disp('Parsed data ok.')
        parsed = true;
    else
        tmp = join(missing_fields(:), ', ');
        disp(['Warning: Parsed data has missing fields: ' tmp{1}] )
    end

    if isfield(parsed_data, 'lanes')
        disp(['Number of detected lanes: ' num2str(length(parsed_data.lanes))])
        tmp = join(parsed_data.lanes(:), ', ');
        disp(['Lanes: ' tmp{1}])
    else
        disp(['WARNING: NO lanes detected. Fix gel_info_file'])
        parsed = false;
    end

end

