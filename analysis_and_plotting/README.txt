This directory contains all analysis scripts for this project. MATLAB scripts (.m files) were used to plot brains using the data generated in R scripts (.R) with the corresponding number (e.g., the results of analyses in "2c - anova BF in windows for brain plots.R" are plotted in "m2c_plot_anova_bfs_on_brains_in_time_bins.m").

"0a - define colors, characters, and ROIs.R" sets up a directory with standard parameters for analyses in this project, including the ROIs, their associated colors, the characters in the stimuli, etc.
"0b - plots ECoG - regions - data warped.R" generates plots in manuscritp Fig. 1b: mean (across trials) of neural activity, per task, per ROI.
"2a - bayesian anovas.R" gets Bayes Factors from Bayesian ANOVAs estimating evidence for word-specific information in each electrode over time (generates data for bottom panel of Fig. 1c).
"2b - anova BF time series by ROI.R" generates plot in bottom panel of Fig. 1c using results of "2a....R"
"2c - anova BF in windows for brain plots.R" generates data for manuscript Fig. 1d: ANOVA Bayes Factors by electrode in time bins.
"m2c_plot_anova_bfs_on_brains_in_time_bins.m" generates brain plots in Fig. 1d using results of "2c....R"
"3c - plot single elec word ECoG high anova samples.R" generates plot in top panel of Fig. 1c: mean neural activity per word.
"7c - word-specific brain plots.R" generates data for Fig. S2: electrodes that are selective for a single word
"m7c_plot_word_specific_elecs.m" uses results of "7c....R" to plot Fig. S2
"8a - cross-validated PN classification - shuffle test  - parallel outer loops.R" trains and tests classifiers on picture naming (PN) data, both for "real" data and 1000 iterations of data with randomly shuffled labels -- the chance distribution.
"8b - get stats from permutations.R" uses the results of "8a....R" to get statistics by comparing the results of "real" data analysis to the distribution of "shuffled" data analysis results.
"8c - plot picture naming cross-validated prediction time series.R" uses the output of "8b....R" to generate plots in Fig. 2b: prediciton accuracy time series for each classifier.
"8f - stack prediction time series.R" stacks the prediction time series plotted in "8c....R" to generate Fig. 2d.
"8g - mean and max time series by ROI.R" gets summary statistics (mean/max) over time of prediction accuracies, and generates the plot in Fig. 2c.
