function folding_screen_enter_lanes(gelData, profileData)
f = figure;

imagesc(gelData.images{1}, [0 3.*std(gelData.images{1}(:))]), axis image, colormap gray, hold on
areas = [profileData.lanePositions(:,1)   profileData.lanePositions(:,3) profileData.lanePositions(:,2)-profileData.lanePositions(:,1)  profileData.lanePositions(:,4)-profileData.lanePositions(:,3)];
for i=1:length(profileData.profiles)
    rectangle('Position', areas(i,:), 'EdgeColor', 'r', 'Linewidth', 1);
    text(areas(i,1)+areas(i,3)/2, areas(i,2) , num2str(i), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 8)
end
set(gca, 'XTickLabel', [], 'YTickLabel', [])


for i=1:length(profileData.profiles)
    c{i} = uicontrol(f,'Style','popupmenu', 'String', 'Bla' );
    c{i}.Position = [0 75+i*25 120 20];
    c{i}.String = {'RM1','RM1_diluted','RM2','M5','M10'};
    c{i}.Callback = @selection;
end


%s=char('TEST')
%set(c{1}.popupmenu1,'string',s)


    function selection(src,event)
        val = src.Value
        str = src.String
        str{val};
        disp(['Selection: ' str{val}]);
    end
disp('bla')
end

%%