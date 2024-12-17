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
band_loop = bands(1)

% Path
main_path = strcat('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/track words during sentences - other versions/gross ROIs with temporal lobe splits along rostral-caudal axis/output/', band_loop, '/');
read_path = strcat(main_path, "data/7c - word-specific brain plots/");

%% Colors
colors_lightest = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_lightest.csv');
colors_lighter = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_lighter.csv');
colors_bright = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_bright.csv');
colors_darker = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_darker.csv');
colors_darkest = readtable('/Users/am4611/Dropbox/Research/ChickenSyntax/analysis/R/color palettes/output/theme_white/rainbow_darkest.csv');
character_colors_t = readtable(strcat(main_path, "data/0a - definitions/character colors/character_colors.csv"));

% Load data and set up colors
word_elecs = struct();
hedges_g = struct();
regions = struct();
locs = struct();
character_colors = struct();
words = string(character_colors_t.word);
for word_loop_n = 1:length(words)
    word_loop = words(word_loop_n);

    % Load word data
    word_elecs.(word_loop) = readtable(strcat(read_path, "hedges g values/elecs selective for just ",word_loop," - hedges g values.csv"));

    % Subset to just elecs in plottable regions
    good_region_indices = find(ismember(string(table2array(word_elecs.(word_loop)(:,"region_clinical"))), plottable_regions));
    word_elecs.(word_loop) = word_elecs.(word_loop)(good_region_indices,:);

    % Remove depth elecs
    non_depth_elec_indices = find(~contains(string(table2array(word_elecs.(word_loop)(:,"elec"))), "_D"))
    word_elecs.(word_loop) = word_elecs.(word_loop)(non_depth_elec_indices ,:);

    % Separate gs, regions, and MNI coordinates
    hedges_g.(word_loop) = table2array(word_elecs.(word_loop)(:,"hedges_g"));
    regions.(word_loop) = string(table2array(word_elecs.(word_loop)(:,"region_clinical")));
    locs.(word_loop) = table2array(word_elecs.(word_loop)(:,["MNI_x","MNI_y","MNI_z"]));

    % Store word color
    character_colors.(word_loop) = table2array(character_colors_t(word_loop_n, ["red","green","blue"])) / 320;
end

% Plot range
min_g = .3;
max_g = .8;


%% Plot

%
% Plot electrode sizes
%

% Save dir
save_cb_dir = fullfile(main_path,...
    "figures",...
    "7c - word-specific brain plots",...
    "color bar - .3 to .8");
if ~exist(save_cb_dir, 'dir')
    mkdir(save_cb_dir);
end

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
linearly_spaced_values = linspace(min_g, max_g, 5);
% Rake radii
elec_size_min = 2.7;
elec_size_max = 5;
fake_radii = ((linearly_spaced_values - min_g) / (max_g - min_g)) * (elec_size_max - elec_size_min) + elec_size_min ;

%     % Get fake alphas
%     value = linearly_spaced_values;
%     c_range = [min(value,[],'all'), w_max];
%     % clip the out of range values
%     value(value < c_range(1)) = c_range(1);
%     value(value > c_range(2)) = c_range(2);
%     % normalize the values
%     rgb_inds = round(((value-c_range(1))/(c_range(2)-c_range(1)))*255)+1;
%     rgb_inds(isnan(rgb_inds)) = 1;
%     % set colors
%     fake_alphas = alpha_map(rgb_inds);

% Plot point sizes
p=VT.PlotElecOnBrain(fake_locations, ...
    'ElecColor', [0,0,0],...
    ...%'ElecAlpha', fake_alphas,...
    'BrainColor', [1,1,1],...
    'cmap',viridis,...
    'clim',[min_g, max_g],...
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
print(gcf, fullfile(save_cb_dir, strcat("elec sizes - .3 to .8.pdf")), '-dpdf', '-bestfit', '-r1000');
close;


%
% Plot brains
%

% Loop through words
for word_loop_n = 1:length(words)
    % word_loop_n = 1
    word_loop = words(word_loop_n)


    %     % Get alphas
    %         value = weights.(cluster_loop);
    %         c_range = [min(value,[],'all'), w_max];
    %         % clip the out of range values
    %         value(value<c_range(1)) = c_range(1);
    %         value(value>c_range(2)) = c_range(2);
    %         % normalize the values
    %         rgb_inds = round(((value-c_range(1))/(c_range(2)-c_range(1)))*255)+1;
    %         rgb_inds(isnan(rgb_inds)) = 1;
    %         % set colors
    %         elec_alphas = alpha_map(rgb_inds);

    % Current weights
    current_regions = regions.(word_loop);
    current_locs = locs.(word_loop);
    current_gs = hedges_g.(word_loop);

    % Cap gs
    current_gs(current_gs > max_g) = max_g;

    % Elec radii:
    elec_radii = ((current_gs - min_g) / (max_g - min_g)) * (elec_size_max - elec_size_min) + elec_size_min ;

    % Project electrodes to nearest pial surface within region
    [elec_locs_proj, ~] = project_elecs_pial(current_locs, current_regions, BrainFile, AnnotFile);

    
    for i = 1:size(current_locs,1)
        project_elecs_pial(current_locs(i,:), current_regions(i), BrainFile, AnnotFile)
        disp(i)
    end
    elec_locs = current_locs(i,:)
    elec_regions = current_regions(i)


    % Initialize brain
    VT = visualtools('Subj', 'MNI',...
        'HS', 'lh',...
        'flag_UseAnnots', false,...
        'BrainFile', BrainFile,...
        'AnnotFile', AnnotFile);

    % Plot elec weights
    p=VT.PlotElecOnBrain(elec_locs_proj, ...
        'ElecColor', repmat(character_colors.(word_loop),[length(current_gs),1]),...
        ...%'ElecAlpha', elec_alphas,...
        'BrainColor', [1,1,1],...%[.6,.6,.6],...
        ...%'cmap',cmap,...
        'clim',[0, max_g],...
        'radius',elec_radii);
    p.AmbientStrength =0.5; %def 0.3
    p.DiffuseStrength =0.6; % def 0.8
    colorbar('off')

    % Set aspect ratio of PDF
    ax = gca; % Get current axes
    aspectRatio = diff(ax.ZLim) / diff(ax.YLim);
    figureWidth = 5; % Define the desired width of the figure in inches
    figureHeight = figureWidth * aspectRatio; % Calculate height to maintain AR
    set(gcf, 'PaperUnits', 'inches'); % Set figure size
    set(gcf, 'PaperSize', [figureWidth figureHeight]);

    % Set up save directory
    save_fig_dir = fullfile(main_path, ...
        "figures",...
        "7c - word-specific brain plots",...
        "brains",...
        "hedges gs - .3 to .8");
    if ~exist(save_fig_dir, 'dir')
        mkdir(save_fig_dir);
    end

    % Save
    print(gcf, fullfile(save_fig_dir, strcat("hedges g - ",word_loop," - .3 to .8.pdf")), '-dpdf', '-bestfit', '-r1000');
    close;
end % word_loop

disp("Script successfully completed.")





