# crisprseeker v1.0.0

What is Compare 2 Sequences?
============================

Compare 2 Sequences generates all possible guide RNAs (gRNAs) for two input sequences, or two sets of sequences and generates scores for potential off-targets in the other sequence.


What is Off Target Analysis?
============================

Off target analysis is the design of target-specific guide RNAs (gRNAs) for CRISPR-Cas9 system by automatically calling findgRNAs, filtergRNAs, searchHits, buildFeatureVectorForScoring, getOfftargetScore, filterOfftarget, calculating gRNA cleavage efficiency and generate reports.

Using the Interface
============================

The user interface on CRISPRSeek allows control over offTargetAnalysis and compare2Sequences by allowing you to change preset defaults depending on your data files. 

First, choose which analysis you would like to do, either off target analysis, or compare 2 sequences. Afterwards, upload your files using the file uploader. If you do not have any files to upload, the default demo files will be used in the analysis instead.

After you have uploaded files, type in the name of the output directory you would like to output the analysis results to. The default output directory is "~/CRISPROut", and will be the output directory unless changed.

Next, change default settings accordingly in the main panel and in the advanced settings panel if needed. Click "Submit" to proceed with analysis, or "Reset all fields to defaults" to reset all the settings you changed.


