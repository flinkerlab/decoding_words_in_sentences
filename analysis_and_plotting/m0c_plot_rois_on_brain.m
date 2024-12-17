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
        colors.(rois(roi_loop_n)) = roi_rgb(roi_loop_n,:);
    end
    
    for brain_color_loop = ["black","white"]
        % brain_color_loop = "white"
    
        if brain_color_loop == "black"
            brain_color = [0,0,0]
        elseif brain_color_loop == "white"
            brain_color = [1,1,1]
        end
    
    %         area_names = {'parsopercularis','parsorbitalis','parstriangularis',... % ifg
%             'caudalmiddlefrontal', 'rostralmiddlefrontal',... % mfg
%             'postcentral','precentral',... % smc
%             'cSTG','mSTG','rSTG',... % stg
%             'cMTG','mMTG','rMTG',... % mtg
%             'supramarginal','inferiorparietal'}; % ipl
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
        
        p=VT.PlotRoIonSurface(area_names, ...
            area_colors / 256, ...
            'BrainColor', brain_color);
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
        save_fig_dir = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/0c - rois on brains/brains/brain legends/');
        if ~exist(save_fig_dir, 'dir')
            mkdir(save_fig_dir);
        end
        
        % Save
        print(gcf, fullfile(save_fig_dir, strcat("rois_cubicl_",brain_color_loop," - ",tl_split_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
        close;


        % %
        % % Plot one ROI at a time
        % %

        for roi_loop_n = 1:length(rois)
            % roi_loop_n = 1
            current_roi = rois(roi_loop_n);
            
            current_areas = roi_components.(current_roi);
            area_colors = repmat(colors.(current_roi), [roi_n_components.(current_roi), 1]);
            
            p=VT.PlotRoIonSurface(current_areas, ...
            area_colors / 255, ...
            'BrainColor', brain_color);
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
            save_individual_fig_dir = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/high_gamma/figures/0c - rois on brains/brains/individual ROIs/');
            if ~exist(save_individual_fig_dir, 'dir')
                mkdir(save_individual_fig_dir);
            end
            
            % Save
            print(gcf, fullfile(save_individual_fig_dir, strcat(current_roi,"_rois_cubicl_",brain_color_loop," - ",tl_split_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
            close;

        end % roi_loop_n

    
    end % brain color
end % tl_split_loop



disp("Script successfully completed.")





