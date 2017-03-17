%%
close all, clear all, clc
addpath('/Users/jonasfunke/MATLAB/MY_MATLAB_TOOLS/dcd_analysis/matdcd-1.0')

fname = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v6_spectrometer/twist_screen/ba_01/output/ba_01-1.dcd';
h = read_dcdheader(fname);

xyz = readdcd(fname, 1:h.N);

%%
% read exb file
exb_fname = '/Users/jonasfunke/Dropbox/FRET_STAGE/Designs/FS-v6_spectrometer/twist_screen/ba_01/ba_01.exb';
fileID = fopen(exb_fname,'r');
C = textscan(fileID,'%s %u %u %f %f ');
fclose(fileID);

extrabond_id = [C{2}, C{3}];
extrabond = [C{4}, C{5}];


%%
size(xyz)

frame = 200;





figure(1), clf
scatter3(xyz(frame,1:3:end), xyz(frame,2:3:end), xyz(frame,3:3:end), 0.3, 'filled'), hold on
%set(hs,'MarkerFaceColor','c');
%alpha(hs,.5);

axis image

%figure(2), clf
for i=1:1000%371 %size(extrabonds,1)
    tmp = [...
        xyz(frame,(extrabond_id(i,1)-45)*3+1), xyz(frame,(extrabond_id(i,1))*3+2), xyz(frame,(extrabond_id(i,1))*3+3); ...
        xyz(frame,(extrabond_id(i,2))*3+1), xyz(frame,(extrabond_id(i,2))*3+2), xyz(frame,(extrabond_id(i,2))*3+3)];
    plot3(tmp(:,1), tmp(:,2), tmp(:,3), 'r-'), hold on
end
    
axis image





