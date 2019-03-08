%%
close all, clear all, clc

filepath = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/JF_Plate_v3/gel_info_simple.txt';
%% read line by line
fileID = fopen(filepath);
tmp = textscan(fileID,'%s','CommentStyle','#', 'Delimiter', '\n');
fclose(fileID);
lines = tmp{1};
disp('File loaded')
out.filename = filepath;
for i=1:length(lines)
    if ~isempty(lines{i}) % if line containes text, e.i. is not empty
        seg = split(lines{i}, {'=', ':'}); % split at = or :
        index_comment = strfind(seg{2}, '#');
        if ~isempty(index_comment)
            seg{2} = seg{2}(1:index_comment(1)-1); % remove all characters afer comment 
        end
        seg{2} = strtrim(seg{2});
        disp(seg{2})
        % user
        if strcmpi(strrep(seg{1}, ' ', ''), 'user')
            out.user = strtrim(seg{2});
        end
        % project
        if strcmpi(strrep(seg{1}, ' ', ''), 'project')
            out.project = strtrim(seg{2});
        end
        % Design_name
        if strcmpi(strrep(seg{1}, ' ', ''), 'design_name')
            out.design_name = strtrim(seg{2});
        end
        % Date
        if strcmpi(strrep(seg{1}, ' ', ''), 'date')
            out.date = strtrim(seg{2});
        end
        % Scaffold_type
        if strcmpi(strrep(seg{1}, ' ', ''), 'scaffold_type')
            out.scaffold_type = strtrim(seg{2});
        end
        if strcmpi(strrep(seg{1}, ' ', ''), 'lattice_type')
            out.lattice_type = strtrim(seg{2});
        end
        if strcmpi(strrep(seg{1}, ' ', ''), 'scaffold_concentration')
            out.scaffold_concentration = str2num(strtrim(seg{2}));
        end
        if strcmpi(strrep(seg{1}, ' ', ''), 'staple_concentration')
            out.staple_concentration = str2num(strtrim(seg{2}));
        end
        if strcmpi(strrep(seg{1}, ' ', ''), 'comment')
            out.comment = strtrim(seg{2});
        end
        
        % Temperature screen
    end
end

out

%% check if all important objects are there
out








