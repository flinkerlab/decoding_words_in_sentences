%% Plot all electrodes from all patients on brain -- project to closest pial surface within anatomical region
clc; clear; close all;

%% Read stuff in
colors = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_bright.csv');
colors_dark = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow.csv');
elecs = readtable("/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/brain plots/ecog/output/data/elec info/patients - combined/row_labels_and_localizations.csv");

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

% Get locations and regions
elec_labels = elecs{:,'patient_elec'};
elec_locs = elecs{:,{'MNI_x','MNI_y','MNI_z'}};
elec_regions = elecs{:,'region_clinical'};


%% Pick elecs
% Elecs with higher SP HGA than LP or PN
elecs_to_plot.sample_anova = ["NY857_G50"];


% Groups
groups = fieldnames(elecs_to_plot);

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
    % color_loop_n = 3
    current_color_label = elec_color_names{color_loop_n};
    current_color_value = elec_colors.(current_color_label);
    disp(current_color_label)

    for group_loop_n = 1:length(groups)
        % group_loop_n = 1
        group_loop = groups{group_loop_n}
        group_elecs = elecs_to_plot.(group_loop);
        
        for elec_loop_n = 1:length(group_elecs)
            % elec_loop_n = 1
            elec_loop = group_elecs(elec_loop_n);
            
            idx = find(strcmp(elec_labels, elec_loop));
            elec_loc = elec_locs(idx,:);
            elec_region = elec_regions(idx);

            if ismember(elec_region, plottable_regions)
    
                % Project electrodes to nearest pial surface within region
                [elec_loc_proj, ~] = project_elecs_pial(elec_loc, elec_region, BrainFile, AnnotFile);
        
                % Bring it forward slightly to make sure it's not partly
                % obscured by some random gyrus
                elec_loc_proj(1) = elec_loc_proj(1) - 20;
        
                % plot electrodes on the brain
                p=VT.PlotElecOnBrain(elec_loc_proj, ...
                                   'ElecColor', current_color_value,...
                                   'BrainColor', [1,1,1],...
                                   'radius',6);
                p.AmbientStrength = 0.5; % default 0.3
                p.DiffuseStrength = 0.6; % default 0.8
                
                % Set aspect ratio of PDF
                ax = gca; % Get current axes
                aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
                figureWidth = 8; % Define the desired width of the figure in inches
                figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
                set(gcf, 'PaperUnits', 'inches'); % Set figure size
                set(gcf, 'PaperSize', [figureWidth figureHeight]);
                
                % Set up save directory
                save_fig_dir = fullfile("/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/3d - plot single elecs on brain/",group_loop,current_color_label);
                if ~exist(save_fig_dir, 'dir')
                    mkdir(save_fig_dir);
                end
                
                % Save
                print(gcf, fullfile(save_fig_dir, strcat(elec_region," - ",elec_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
                
                close;
            else % if elec is in a plottable region
                disp(strcat(elec_loop, " is not in a plottable region (",elec_region,"). Skipping."))
            end % if elec is in a plottable region
        end % elec_loop
    end % group_loop
end % color_loop

disp("Script completed successfully!")