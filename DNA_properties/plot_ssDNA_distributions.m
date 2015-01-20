%%
close all, clear all, clc

%%
n= 5:5:100;
cc = varycolor(length(n));
n_mean = ssDNA_mean_end_to_end_distance(n, 2.3); % length of leash in bases, persitence length in bases
close all
hold all
for i=1:length(n)
    tmp = ssDNA_end_to_end_distribution(n(i),  2.3);
    plot(tmp(:,1), tmp(:,2), 'Color', cc(i,:))
    set(gca, 'XLim', [0 40])
    vline(n_mean(i), {'Color', cc(i,:)});

    pause
end


%%
