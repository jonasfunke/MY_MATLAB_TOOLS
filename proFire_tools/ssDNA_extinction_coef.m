function [e] = ssDNA_extinction_coef(seq)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% make sequence all captal
disp(['Input sequence:' seq])
seq = upper(seq);
seq = strrep(seq,' ','');
seq_del = seq;
seq_del = strrep(seq_del,'G','');
seq_del = strrep(seq_del,'T','');
seq_del = strrep(seq_del,'C','');
seq_del = strrep(seq_del,'A','');

    if ~isempty(seq_del)
        disp(['ERROR: Non-DNA bases detected. Only use, G,C,T, and A. Non-DNA characters found: ' seq_del])
        e = 0;
    else
        % number of single bases
        N(1) = sum(seq=='A');
        N(2) = sum(seq=='C');
        N(3) = sum(seq=='G');
        N(4) = sum(seq=='T');

        e_bases = [15400 7400 11500 8700];

        %% number of nearest neighbours
        NN(1,1) = length(strfind(seq,'AA')); 
        NN(1,2) = length(strfind(seq,'AC')); 
        NN(1,3) = length(strfind(seq,'AG')); 
        NN(1,4) = length(strfind(seq,'AT')); 

        NN(2,1) = length(strfind(seq,'CA')); 
        NN(2,2) = length(strfind(seq,'CC')); 
        NN(2,3) = length(strfind(seq,'CG')); 
        NN(2,4) = length(strfind(seq,'CT'));

        NN(3,1) = length(strfind(seq,'GA')); 
        NN(3,2) = length(strfind(seq,'GC')); 
        NN(3,3) = length(strfind(seq,'GG')); 
        NN(3,4) = length(strfind(seq,'GT'));

        NN(4,1) = length(strfind(seq,'TA')); 
        NN(4,2) = length(strfind(seq,'TC')); 
        NN(4,3) = length(strfind(seq,'TG')); 
        NN(4,4) = length(strfind(seq,'TT'));

        e_NN = [27400 21200 25000 22800; ...
            21200 14600 18000 15200; ...
            25200 17600 21600 20000; ...
            23400 16200 19000 16800];

        % calculate extinction coef
        e = sum(sum(NN.*e_NN))-sum(N.*e_bases);

    end
end
