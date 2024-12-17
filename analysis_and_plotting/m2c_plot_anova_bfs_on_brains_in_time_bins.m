%% Plot all electrodes from all patients on brain -- project to closest pial surface within anatomical region
clc; clear; close all;

%% VisualTools stuff
cd '/Users/am4611/Dropbox/Research/ChickenSyntax/brain plot code from amir/visualization-tools-main/matlab/'
%AnnotFile='./SampleData/MNI/ch2_template.lh.aparc.split_STG_MTG.annot';
%BrainFile='./SampleData/MNI/ch2_template_lh_pial_120519.mat';
AnnotFile='./SampleData/MNI-FS/FSL_MNI152.lh.aparc.split_STG_MTG.annot';
BrainFile='./SampleData/MNI-FS/FSL_MNI152_lh_pial.mat';
[~, ~, plottable_regions] = visualtools.read_annot(AnnotFile);
plottable_regions = string(plottable_regions.struct_names);


%% Metadata
% Band loops
bands = ["high_gamma", "beta"];
band_loop_n = 1;
band_loop = bands(band_loop_n)

% Path
main_path = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/', band_loop, '/');
path_to_bin_sizes = strcat(main_path, 'data/2c - anova BF in windows for brain plots/');
bin_sizes = string({dir(path_to_bin_sizes).name});
bin_sizes = bin_sizes(~ismember(bin_sizes, {'.','..','.DS_Store'}));

for bin_size_loop_n = 1:length(bin_sizes)
    % bin_size_loop_n = 1

    bin_size_loop = bin_sizes(bin_size_loop_n);

    datasets = ["elec_max_log_BF_in_bins", "elec_mean_log_BF_in_bins"];

    for dataset_loop_n = 1:length(datasets)
        % dataset_loop_n = 1

        dataset_loop = datasets(dataset_loop_n);

        % Load data
        data = readtable(strcat(path_to_bin_sizes,bin_size_loop,'/',dataset_loop,'.csv'));

        % Get rid of regions not in VisualTools (mostly depth elecs and NAs)
        good_region_indices = find(ismember(string(table2array(data(:,"region_clinical"))), plottable_regions));
        data = data(good_region_indices,:);

        % Separate data from metadata
        window_columns = startsWith(data.Properties.VariableNames, 'bin');
        elecs = string(data.elec);
        regions = string(data.region_clinical);
        locs = table2array(data(:,["MNI_x","MNI_y","MNI_z"]));
        data = data(:, window_columns);

        %% Colors
        colors_lightest = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_lightest.csv');
        colors_lighter = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_lighter.csv');
        colors_bright = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_bright.csv');
        colors_darker = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_darker.csv');
        colors_darkest = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_darkest.csv');
        
        
        % Color max val
        lower_cap = 0;
        upper_cap = 8;
        
        %% Loop thru colors
        color_labels = "pink";%["purple", "orange", "green", "pink", "cyan", "blue","red"];
        
        %for color_loop_n = 1:length(color_labels)
        color_loop_n = 1
        color_label = color_labels(color_loop_n);
        color_scale = [
            table2array(colors_darkest(strcmp(color_label, colors_darkest.Var1),2:4));
            %table2array(colors_darker(strcmp(color_label, colors_darker.Var1),2:4));
            %table2array(colors_bright(strcmp(color_label, colors_bright.Var1),2:4));
            table2array(colors_lighter(strcmp(color_label, colors_lighter.Var1),2:4));
            %table2array(colors_lightest(strcmp(color_label, colors_lightest.Var1),2:4));
            1,1,1];
        alpha_scale = [1,1,1,.8,0]';
    
    
        cmap = interp1(linspace(lower_cap, upper_cap, size(color_scale, 1)), color_scale, linspace(upper_cap, 0, 256),'pchip');
        alpha_map = interp1(linspace(lower_cap, upper_cap, size(alpha_scale, 1)), alpha_scale, linspace(upper_cap, 0, 256),'pchip');
    
        %
        % Plot dummy data and create colorbar
        %
    
        figure;
        ax = axes('Visible', 'off');
        fake_data = zeros(10);
        fake_data(1,1) = upper_cap;
        fake_data(1,2) = lower_cap;
        plt = imagesc(ax, fake_data); % Plot dummy data
        colormap(ax, cmap); % Set the colormap you want to use
        c = colorbar;
    
        % Adjust
        plt.Visible = 'on';
        ax.Visible = 'off';
        c.Visible = 'on'; % Make only the colorbar visible
    
        % Optionally set colorbar properties
        set(c, 'Box', 'off');
        set(c, 'TickLength', 0);
    
        % Save colorbar
        save_cb_dir = fullfile(main_path,"figures","2c - anova BF in windows for brain plots","brains","brain color bar");
        if ~exist(save_cb_dir, 'dir')
            mkdir(save_cb_dir);
        end
        print(gcf, fullfile(save_cb_dir, strcat("color bar.pdf")), '-dpdf', '-bestfit', '-r1000');
        close;
    
        %
        % Plot electrode sizes
        %
    
        VT = visualtools('Subj', 'MNI',...
            'HS', 'lh',...
            'flag_UseAnnots', false,...
            'BrainFile', BrainFile,...
            'AnnotFile', AnnotFile);
    
        % Plot elec weights
        fake_locations = [0,0,-40; ...
            0,0,-20; ...
            0,0,0; ...
            0,0,20; ...
            0,0,40];
        weight_min_max = [lower_cap, upper_cap];
        linearly_spaced_values = linspace(weight_min_max(1), weight_min_max(2), 5);
        % Fake radii
        elec_size_min = 1;
        elec_size_max = 2.7;
        fake_radii = ((linearly_spaced_values - lower_cap) / (upper_cap - lower_cap)) * (elec_size_max - elec_size_min) + elec_size_min ;
        
        % Get fake alphas
        values = linearly_spaced_values;
        c_range = [lower_cap, upper_cap];
        % clip the out of range values
        values(values < c_range(1)) = c_range(1);
        values(values > c_range(2)) = c_range(2);
        % normalize the values
        rgb_inds = round(((values-c_range(1))/(c_range(2)-c_range(1)))*255)+1;
        rgb_inds(isnan(rgb_inds)) = 1;
        % set colors
        fake_alphas = alpha_map(rgb_inds);
    
        % Plot point sizes
        p=VT.PlotElecOnBrain(fake_locations, ...
            'ElecColor', linearly_spaced_values',...
            'ElecAlpha', fake_alphas,...
            'BrainColor', [1,1,1],...
            'cmap',cmap,...
            'clim',[lower_cap, upper_cap],...
            'radius',fake_radii, ...
            'FaceAlpha',0,...
            'flag_JustElecs',false);

        % Set aspect ratio of PDF
        ax = gca; % Get current axes
        aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
        figureWidth = 5; % Define the desired width of the figure in inches
        figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
        set(gcf, 'PaperUnits', 'inches'); % Set figure size
        set(gcf, 'PaperSize', [figureWidth figureHeight]);
    
        % Save elec sizes
        print(gcf, fullfile(save_cb_dir, strcat("elec sizes.pdf")), '-dpdf', '-bestfit', '-r1000');
        close;
    
    
        %% Loop thru windows
        windows = string(data.Properties.VariableNames);
        for(window_loop_n = 1:length(windows))
            % window_loop_n = 1
            window_loop = windows(window_loop_n)
    
            % Get alphas
            values = data.(window_loop);
            c_range = [lower_cap, upper_cap];
            % clip the out of range values
            values(values<c_range(1)) = c_range(1);
            values(values>c_range(2)) = c_range(2);
            % normalize the values
            rgb_inds = round(((values-c_range(1))/(c_range(2)-c_range(1)))*255)+1;
            rgb_inds(isnan(rgb_inds)) = 1;
            % set colors
            elec_alphas = alpha_map(rgb_inds);

            % Elec radii:
            elec_radii = ((values - lower_cap) / (upper_cap - lower_cap)) * (elec_size_max - elec_size_min) + elec_size_min ;

            % Project electrodes to nearest pial surface within region
            [elec_locs_proj, ~] = project_elecs_pial(locs, regions, BrainFile, AnnotFile);

            % Initialize brain
            VT = visualtools('Subj', 'MNI',...
                'HS', 'lh',...
                'flag_UseAnnots', false,...
                'BrainFile', BrainFile,...
                'AnnotFile', AnnotFile);

            % Plot only elecs with BFs > 3
            bf_over_3_indices = values > log10(3);

%             % Plot coverage
%             p=VT.PlotElecOnBrain(elec_locs_proj, ...
%                 'ElecColor', [.5, .5, .5],...
%                 ...%'ElecAlpha', elec_alphas,...
%                 'BrainColor', [1,1,1],...
%                 'cmap',cmap,...
%                 'clim',[lower_cap, upper_cap],...
%                 'radius',.7);
%             p.AmbientStrength =0.5; %def 0.3
%             p.DiffuseStrength =0.6; % def 0.8
            
            % Plot elec weights
            p=VT.PlotElecOnBrain(elec_locs_proj(bf_over_3_indices,:), ...
                'ElecColor', values(bf_over_3_indices),...
                ...%'ElecAlpha', elec_alphas,...
                'BrainColor', [1,1,1],...
                'cmap',cmap,...
                'clim',[lower_cap, upper_cap],...
                'radius',elec_radii(bf_over_3_indices), ...
                'flag_JustElecs',false);
            p.AmbientStrength =0.5; % def 0.3
            p.DiffuseStrength =0.6; % def 0.8
            
            % Turn off color bar
            colorbar('off')

            % Set aspect ratio of PDF
            ax = gca; % Get current axes
            aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
            % figureWidth = 5; % DEFINED ABOVE for "elec sizes" plot (must
            % be identical)
            figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
            set(gcf, 'PaperUnits', 'inches'); % Set figure size
            set(gcf, 'PaperSize', [figureWidth figureHeight]);

            % Set up save directory
            save_fig_dir = fullfile(main_path, ...
                "figures",...
                "2c - anova BF in windows for brain plots",...
                "brains",...
                dataset_loop,...
                bin_size_loop);
            if ~exist(save_fig_dir, 'dir')
                mkdir(save_fig_dir);
            end

            % Save
            print(gcf, fullfile(save_fig_dir, strcat(window_loop," - ",dataset_loop,".pdf")), '-dpdf', '-bestfit', '-r1000');
            close;
        end % window_loop_n
        % end % color_loop_n
    end % dataset_loop_n
end % bin_size_loop_n

disp("Script successfully completed.")





