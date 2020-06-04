%%


fid=fopen('200601_1248_Cy3_test_2_differentEx_sameEm.csv');
tline = fgetl(fid);
while ischar(tline)
    %tlines{end+1,1} = tline;
    tline = fgetl(fid);
    %disp(tline)
    if startsWith(tline, ' Number of cycles:')
        tmp = textscan(tline(19:end),'%d');
        N_cycles = tmp{1};
    end
    
    if startsWith(tline, 'Time')
        display(['Time detected...' ])
        time = fgetl(fid) % read time

        tline = fgetl(fid);
        while ischar(tline) && ~startsWith(tline, 'Chromatic')
            pos = tline(1:3);
            tmp = textscan(tline(5:end),'%d');
            %disp(tline); %write to file
            tline = fgetl(fid);
        end

    end
end
fclose(fid);