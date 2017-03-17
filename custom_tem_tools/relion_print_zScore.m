%%

close all, clear all, clc

%%
fileID = fopen('/Users/jonasfunke/Dropbox/Temporary/2016-08-22/particles_autopick_all_ctf.txt');
C = textscan(fileID,'%s %.6f %.6f %s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %u %.6f','MultipleDelimsAsOne',1, 'HeaderLines',21);
fclose(fileID);

%%
zScores =C{17};

%%
zScores = sort(zScores);

%%

figure(1);
plot(zScores)
%set(gca, 'YLim', [0 1])
%C = textscan(fileID,formatSpec,N,'CommentStyle','##','Delimiter','\t');
