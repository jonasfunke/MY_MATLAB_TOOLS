function parsed_data = parse_gel_info_simple(filepath)

    %% read line by line
    fileID = fopen(filepath);
    tmp = textscan(fileID,'%s','CommentStyle','#', 'Delimiter', '\n');
    fclose(fileID);
    lines = tmp{1};
    parsed_data.filename = filepath;
    for i=1:length(lines)
        if ~isempty(lines{i}) % if line containes text, e.i. is not empty
            seg = split(lines{i}, {'=', ':'}); % split at = or :
            index_comment = strfind(seg{2}, '#');
            if ~isempty(index_comment)
                seg{2} = seg{2}(1:index_comment(1)-1); % remove all characters afer comment 
            end
            seg{2} = strtrim(seg{2});
            %disp(seg{2})
            % user
            if strcmpi(strrep(seg{1}, ' ', ''), 'user')
                parsed_data.user = strtrim(seg{2});
            end
            % project
            if strcmpi(strrep(seg{1}, ' ', ''), 'project')
                parsed_data.project = strtrim(seg{2});
            end
            % Design_name
            if strcmpi(strrep(seg{1}, ' ', ''), 'design_name')
                parsed_data.design_name = strtrim(seg{2});
            end
            % Date
            if strcmpi(strrep(seg{1}, ' ', ''), 'date')
                parsed_data.date = strtrim(seg{2});
            end
            % Scaffold_type
            if strcmpi(strrep(seg{1}, ' ', ''), 'scaffold_type')
                parsed_data.scaffold_type = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'lattice_type')
                parsed_data.lattice_type = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'scaffold_concentration')
                parsed_data.scaffold_concentration = str2num(strtrim(seg{2}));
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'staple_concentration')
                parsed_data.staple_concentration = str2num(strtrim(seg{2}));
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'comment')
                parsed_data.comment = strtrim(seg{2});
            end

            % Lanes
            for l=1:20
                if strcmpi(strrep(seg{1}, ' ', ''), ['Lane_' sprintf('%02i', l)])
                    parsed_data.lanes{l} = strtrim(seg{2});
                end
            end

            % Gel parameters
            if strcmpi(strrep(seg{1}, ' ', ''), 'Gelsize')
                parsed_data.gelsize = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'Agarose_concentration')
                parsed_data.agarose_concentration = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'Staining')
                parsed_data.staining = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'Mg_concentration')
                parsed_data.mg_concentration = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'Voltage')
                parsed_data.voltage = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'Running_time')
                parsed_data.running_time = strtrim(seg{2});
            end
            if strcmpi(strrep(seg{1}, ' ', ''), 'Cooling')
                parsed_data.cooling = strtrim(seg{2});
            end
        end
    end
    disp(['File ' filepath ' parsed.'])
end








