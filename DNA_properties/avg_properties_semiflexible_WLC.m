function [ r_values ] = avg_properties_semiflexible_WLC( l_c, l_p, display )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    dL = l_c/1000;

    %distribution = cell(length(L),1);

        if l_c <= 0
            r_mean = 0;
            r_std = 0;
            r_median = 0;
        else
            R = 0:dL:l_c-dL;
            x = l_p*(1-R./l_c)./l_c;
            end_to_end = (l_p/l_c^2).*frey_function(x);
            end_to_end = end_to_end./sum(end_to_end);
            
           % distribution = [R', end_to_end'];
            
            r_mean = sum(end_to_end.*R); % mean = r1
            r2 = sum(end_to_end.*R.^2); % mean
            r_std = sqrt( r2 - r_mean.^2  );
            
            %calculate median 
            cumulant = cumsum(end_to_end);
            R_sub = R( cumulant < 0.6 & cumulant > 0.4);
            c_sub = cumulant( cumulant < 0.6 & cumulant > 0.4);
            
            r_median = spline(c_sub, R_sub, 0.5);
            
            if display
                %plot(R, end_to_end, R, cumulant.*max(end_to_end), R_sub, c_sub.*max(end_to_end)), hold on
                plot(R, end_to_end./max(end_to_end)), hold on
               % plot(R, cumulant.*max(end_to_end)), hold on
               % vline(r_mean, 'k');
               % vline([r_mean-r_std, r_mean+r_std], 'k--');
               % vline(r_median, 'k:');
            end
            
             r_values = [r_mean, r_std, r_median];
        end

end

