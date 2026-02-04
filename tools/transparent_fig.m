
% Exporting figures with transparent background

base_folder = "C:\Users\Ariosky\OneDrive - CCLAB\Data\Ludovico\Outputs\structural";
figures = dir(fullfile(base_folder,'*.fig'));
for i=1:length(figures)
    figure = figures(i);
    fig = openfig(fullfile(figure.folder,figure.name));
    [~,fig_name,~] = fileparts(figure.name);
    views = {'left','right','front','back','top'};
    for v=1:length(views)
        rotate_fig(views{v});
        export_fig(fig,fullfile(figure.folder,strcat(fig_name,'-',views{v})),'-transparent','-png');
    end
    close(fig);
end




