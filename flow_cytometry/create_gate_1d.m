function [ I_gate ] = create_gate_1d( I_in, channel_name, message )

%  let the user select a gate value
%  Input: I= liste (1d) of intensities
%  Output: I_gate = gate valie
    
    I = real(log10(I_in));
    xlim = [min(I) (max(I))];
    ylim = [0 0];
    
    cur_fig = figure;
    tmp = histogram(I, 'DisplayStyle', 'stairs'); hold on
    if max(tmp.Values)>ylim(2)
        ylim(2) = max(tmp.Values);
    end
    title(message)
    xlabel(channel_name)
    ylabel('Counts')
    set(gca, 'XLim', xlim, 'YLim', ylim)

    [x_selected,~] = ginput(1);
    I_gate = 10^x_selected;
    close(cur_fig)


end

