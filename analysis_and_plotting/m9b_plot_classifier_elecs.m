%% Plot all electrodes from all patients on brain -- project to closest pial surface within anatomical region
clc; clear; close all;

%% Read stuff in
colors = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_bright.csv');
colors_dark = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow.csv');

% Get colors
elec_colors.pink = table2array(colors(strcmp("pink", colors.Var1),2:4));
%elec_colors.blue = table2array(colors(strcmp("blue", colors.Var1),2:4));
elec_colors.black = [0,0,0];
elec_color_names = fieldnames(elec_colors);

% Brain color parameters specific to colors
brain_color.pink = [1,1,1];
brain_color.blue = [1,1,1];
brain_color.black = [1,1,1];

brain_ambient_strength.pink = .5;
brain_ambient_strength.blue = .5;
brain_ambient_strength.black = .5;

brain_diffuse_strength.pink = .6;
brain_diffuse_strength.blue = .6;
brain_diffuse_strength.black = .6;

%% Pick elecs
% Elecs with higher SP HGA than LP or PN
classifiers_to_plot.fig2 = ["IFG_NY863",
    "SMC_NY837",
    "PTL_NY863",
    "MTL_NY765",
    "SMC_NY829"];
classifiers_to_plot.fig3 = ["SMC",
    "ATL"];


% Groups
groups = fieldnames(classifiers_to_plot);

%% Plot!
cd '/Users/am4611/Dropbox/Research/ChickenSyntax/brain plot code from amir/visualization-tools-main/matlab/'
%AnnotFile='./SampleData/MNI/ch2_template.lh.aparc.split_STG_MTG.annot';
%BrainFile='./SampleData/MNI/ch2_template_lh_pial_120519.mat';
AnnotFile='./SampleData/MNI-FS/FSL_MNI152.lh.aparc.split_STG_MTG.annot';
BrainFile='./SampleData/MNI-FS/FSL_MNI152_lh_pial.mat';
[~, ~, plottable_regions] = visualtools.read_annot(AnnotFile);
plottable_regions = string(plottable_regions.struct_names);

% Initialize brain
VT = visualtools('Subj', 'MNI',...
                 'HS', 'lh',...
                 'flag_UseAnnots', false,...
                 'BrainFile', BrainFile,...
                 'AnnotFile', AnnotFile);

%% Plot elec means
for color_loop_n = 1:length(elec_color_names)
    % color_loop_n = 1
    current_color_label = elec_color_names{color_loop_n};
    current_color_value = elec_colors.(current_color_label);
    disp(current_color_label)

    for group_loop_n = 1:length(groups)
        % group_loop_n = 1
        group_loop = groups{group_loop_n};
        group_classifiers = classifiers_to_plot.(group_loop);
        
        for classifier_loop_n = 1:length(group_classifiers)
            % classifier_loop_n = 1
            classifier_loop = group_classifiers(classifier_loop_n);
            
            % Load elecs in this classifier
            classifier_elecs_t = readtable(strcat("/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/data/9b - get classifier elecs for brain plots/classifier elecs dataframes/",classifier_loop,"_elecs.csv"));
            
            % Limit to good rows -- non depth elecs in plottable regions
            good_elec_indices = find(~contains(classifier_elecs_t.elec,"_D"));
            classifier_elecs_t = classifier_elecs_t(good_elec_indices,:);
            good_region_indices = find(ismember(string(table2array(classifier_elecs_t(:,"region_clinical"))), plottable_regions));
            classifier_elecs_t = classifier_elecs_t(good_region_indices,:);
            
            % Get plot data
            locs = table2array(classifier_elecs_t(:,["MNI_x","MNI_y","MNI_z"]));
            regions = string(table2array(classifier_elecs_t(:,"region_clinical")));

            % Project electrodes to nearest pial surface within region
            [elec_loc_proj, ~] = project_elecs_pial(locs, regions, BrainFile, AnnotFile);
    
    
            % plot electrodes on the brain
            p=VT.PlotElecOnBrain(elec_loc_proj, ...
                               'ElecColor', current_color_value,...
                               'BrainColor', [1,1,1],...
                               'radius',4.5);
            p.AmbientStrength = 0.5; % default 0.3
            p.DiffuseStrength = 0.6; % default 0.8
            
            % Set aspect ratio of PDF
            ax = gca; % Get current axes
            aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
            figureWidth = 5; % Define the desired width of the figure in inches
            figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
            set(gcf, 'PaperUnits', 'inches'); % Set figure size
            set(gcf, 'PaperSize', [figureWidth figureHeight]);
            
            % Set up save directory
            save_fig_dir = fullfile("/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/9b - plot classifier elecs on brain/",group_loop,current_color_label);
            if ~exist(save_fig_dir, 'dir')
                mkdir(save_fig_dir);
            end
            
            % Save
            print(gcf, fullfile(save_fig_dir, strcat(classifier_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
            
            close;
        
        end % elec_loop
    end % group_loop
end % color_loop

disp("Script completed successfully!")