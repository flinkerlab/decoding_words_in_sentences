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
"9a - train classifiers - get variable importances - no cross-validation - logistic regression.R" retrains all classifiers using 100% of the picture naming data (rather than 90% as in previous analyses to facilitate 10-fold cross validation).
"10a - decode sentence nouns - shuffle test.R" uses the classifiers trained in "9a....R" to predict word identity during the production of lists, active sentences, and passive sentences. Performed once for "real" data and 1000 times for shuffled data, as in "8b....R".
"10b_vE - get sentence stats from permutations - empirical quantile - vD but p<.05 for 100ms or p<.01 for 50ms.R" gets stats from the permutation distributions generated in "10a....R".
"12b - plot sentence prediction time series stack.R" generates stack plots in Figs. 3b, 4b, and 4d (corresponding to time series plots generated by "18c_vA....R").
"15a - barplots - proportion classifiers sig.R" counts the number of "congruent" and "incongruent" noun detections and generates barplots in Fig. 4e.
"18c_vA - plot sentence and list prediction time series.R" generates the time series plots of prediction accuracy in Figs. 3a, 4a, and 4c (corresponding to the stacked time series plots generated in "12b....R").
"21a - plot significant prediction densities by ROI.R" generates density plots in Fig. S4.
