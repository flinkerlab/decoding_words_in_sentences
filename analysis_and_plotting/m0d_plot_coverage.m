%% Plot all electrodes from all patients on brain -- project to closest pial surface within anatomical region
clc; clear; close all;

% VisualTools stuff
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
                 'flag_UseAnnots', true,...
                 'BrainFile', BrainFile,...
                 'AnnotFile', AnnotFile);


roi_components.IFG = ["parsopercularis","parsorbitalis","parstriangularis"];
roi_components.SMC = ["precentral","postcentral"];
roi_components.ATL = ["temporalpole","rMTG","rSTG"];
roi_components.MTL = ["mMTG","mSTG"];
roi_components.PTL = ["cMTG","cSTG"];
roi_components.IPL = ["supramarginal","inferiorparietal"];
roi_components.STG = ["cSTG","mSTG","rSTG"]; % "ctx-lh-superiortemporal"
roi_components.MTG = ["cMTG","mMTG","rMTG"]; % "ctx-lh-middletemporal"
roi_components.MFG = ["caudalmiddlefrontal","rostralmiddlefrontal"];

% Load elec locations
elec_info = readtable("/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/data/0a - definitions/elec info/elec info.csv");

% Get elecs in plottable regions
keep_region_indices = find(ismember(string(elec_info.region_clinical_for_plots), plottable_regions) & (~ startsWith(string(elec_info.elec), "D")));

% Subset and get data
elec_info = elec_info(keep_region_indices,:);
locs = table2array(elec_info(:,["MNI_x","MNI_y","MNI_z"]));
regions = string(elec_info.region_clinical_for_plots);

% Project electrodes to nearest pial surface within region
[elec_locs_proj, ~] = project_elecs_pial(locs, regions, BrainFile, AnnotFile);


% %
% % Plot just brain
% % 

tl_splits = ["dorsal.ventral", "rostral.caudal"];

for tl_split_loop_n = 1:2
    % tl_split_loop_n = 1
    tl_split_loop = tl_splits(tl_split_loop_n);
    
    roi_table = readtable(strcat("/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/data/0a - definitions/roi metadata/roi colors with temporal lobe split along ",tl_split_loop," axis.csv"))
    rois = string(roi_table.roi);
    roi_rgb = table2array(roi_table(:,["red","green","blue"]));

    colors = struct();
    for roi_loop_n = 1:length(rois)
        colors.(rois(roi_loop_n)) = roi_rgb(roi_loop_n,:) / 255;
    end
    colors.None = [.5, .5, .5];
    color_names = string(fieldnames(colors));

    % Make a matrix of colors for elecs based on region
    elec_roi_colors = nan(length(regions), 3);
    for elec_loop = 1:length(regions)
        current_roi = "None";
        for i = 1:length(rois)
            if ismember(regions(elec_loop), roi_components.(rois(i)));
                current_roi = color_names(i);
                break;
            end
        end
        elec_roi_colors(elec_loop,:) = colors.(current_roi);
    end % elec_loop

    
    for brain_color_loop = ["black","white"]
        % brain_color_loop = "white"
    
        if brain_color_loop == "black"
            brain_color = [0,0,0]
            elec_color = [.999,.999,.999]
        elseif brain_color_loop == "white"
            brain_color = [1,1,1]
            elec_color = [0,0,0]
        end
    
        area_names = [];
        for roi_loop_n = 1:length(rois)
            current_roi = rois(roi_loop_n);
            area_names = [area_names, roi_components.(current_roi )];
            roi_n_components.(current_roi ) = length(roi_components.(current_roi ));
        end
        area_names = area_names';
        
        area_colors = [];
        for roi_loop_n = 1:length(rois)
            current_roi = rois(roi_loop_n);
            area_colors = [area_colors; repmat(colors.(current_roi), [roi_n_components.(current_roi), 1])];
        end
        
        
        %
        % Elecs on colored regions
        %

        VT = visualtools('Subj', 'MNI',...
                 'HS', 'lh',...
                 'flag_UseAnnots', true,...
                 'BrainFile', BrainFile,...
                 'AnnotFile', AnnotFile);

        p=VT.PlotRoIonSurface(area_names, ...
            area_colors, ...
            'BrainColor', brain_color);
        p.AmbientStrength =0.5; %def 0.3
        p.DiffuseStrength =0.6; % def 0.8

        % Add elecs
        p=VT.PlotElecOnBrain(elec_locs_proj, ...
                'ElecColor', elec_color,...
                ...%'BrainColor', [1,1,1],...
                'radius',1, ...
                'flag_JustElecs',true);
        
        % Set aspect ratio of PDF
        ax = gca; % Get current axes
        aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
        figureWidth = 5; % Define the desired width of the figure in inches
        figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
        set(gcf, 'PaperUnits', 'inches'); % Set figure size
        set(gcf, 'PaperSize', [figureWidth figureHeight]);
        
        % Set up save directory
        save_fig_dir = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/0d - elec coverage and brain regions/brains/elecs on colored regions/');
        if ~exist(save_fig_dir, 'dir')
            mkdir(save_fig_dir);
        end
        
        % Save
        print(gcf, fullfile(save_fig_dir, strcat("rois_cubicl_",brain_color_loop," - ",tl_split_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
        close;



        %
        % Elecs on plain brain
        %

        % Only do this once (should be same on both tl_split_loops)
        if tl_split_loop_n == 2

            VT = visualtools('Subj', 'MNI',...
                'HS', 'lh',...
                'flag_UseAnnots', false,...
                'BrainFile', BrainFile,...
                'AnnotFile', AnnotFile);
    
            p=VT.PlotElecOnBrain(elec_locs_proj, ...
                    'ElecColor', elec_color,...
                    'BrainColor', brain_color,...
                    ...%'cmap',cmap,...
                    ...%'clim',[lower_cap, upper_cap],...
                    'radius',1);
            p.AmbientStrength =0.5; %def 0.3
            p.DiffuseStrength =0.6; % def 0.8
    
            
            % Set aspect ratio of PDF
            ax = gca; % Get current axes
            aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
            figureWidth = 5; % Define the desired width of the figure in inches
            figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
            set(gcf, 'PaperUnits', 'inches'); % Set figure size
            set(gcf, 'PaperSize', [figureWidth figureHeight]);
            
            % Set up save directory
            save_fig_dir = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/0d - elec coverage and brain regions/brains/elecs on plain brain/');
            if ~exist(save_fig_dir, 'dir')
                mkdir(save_fig_dir);
            end
            
            % Save
            print(gcf, fullfile(save_fig_dir, strcat("rois_cubicl_",brain_color_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
            close;

        end %if tl_split_loop_n == 2


        %
        % Colored elecs on plain brain
        %

        VT = visualtools('Subj', 'MNI',...
            'HS', 'lh',...
            'flag_UseAnnots', false,...
            'BrainFile', BrainFile,...
            'AnnotFile', AnnotFile);

        p=VT.PlotElecOnBrain(elec_locs_proj, ...
                'ElecColor', elec_roi_colors,...
                'BrainColor', brain_color,...
                ...%'cmap',cmap,...
                ...%'clim',[lower_cap, upper_cap],...
                'radius',1.4);
        p.AmbientStrength =0.5; %def 0.3
        p.DiffuseStrength =0.6; % def 0.8

        
        % Set aspect ratio of PDF
        ax = gca; % Get current axes
        aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
        figureWidth = 5; % Define the desired width of the figure in inches
        figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
        set(gcf, 'PaperUnits', 'inches'); % Set figure size
        set(gcf, 'PaperSize', [figureWidth figureHeight]);
        
        % Set up save directory
        save_fig_dir = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/0d - elec coverage and brain regions/brains/colored elecs on plain brain/');
        if ~exist(save_fig_dir, 'dir')
            mkdir(save_fig_dir);
        end
        
        % Save
        print(gcf, fullfile(save_fig_dir, strcat("rois_cubicl_",brain_color_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
        close;
        
    
    end % brain color
end % tl_split_loop



disp("Script successfully completed.")





