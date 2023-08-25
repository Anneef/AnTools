<!DOCTYPE html>
<html>
<body>
  
# Overview
This Igor Pro was developed by Andreas Neef for Igor Pro 7 to 9. It falls into three categories. *AnTools_extended.ipf* was developed to automate general wave- and string-based tasks as well as basic data input tasks and graphing tasks.
The collection inside *Dynamic Gain code* contains some general tools to analyze cellular electrophysiology, specifically spike detection and analysis, and the analysis of spikes fired in response to conventional, square pulse stimulation (frequency-current curves etc.). *Dynamic Time Warping* is the translation of Matlab code from David Schultz (see below).
  
More importantly, the Dynamic Gain code comprises all code to compute the dynamic gain function of neurons that fire in response to continuously fluctuating, <i>in vivo</i>-like stimuli. This analysis reveals how populations encode information under realistic stimuli, aspects of which cannot be obtained with conventional stimuli ( [Lazarov et al. 2018 Sci.Advances](https://doi.org/10.1126/sciadv.aau8621), [Revah et al. 2019 J Neurosci.](https://doi.org/10.1523/JNEUROSCI.3147-18.2019) and [Merino et al. 2021 PNAS](https://doi.org/10.1073/pnas.2114549118)).

## AnTools
ANTools_extended.ipf contains useful tools to speed up common steps in general data analysis and data display.

### *i) Wrapper functions for repeated application of a list of commands for each item.*
  
  Inside the commands,
   a number of wildcards (~ §) stand in for the wave name, the index of the current wave in the list of waves etc. Those wildcards are elaborated in the beginning of the code for the respective functions.
   The items that are looped over can be 
  - a list of waves provided directly - with extra options
    
    `Xeqt4List("WaveName1; Wavename2;","ExcludeThisName;ExcludeThatNameToo;",Boolean-UseFullPathInCommand,"Cmd1 ~; Cmd2 ~;")`
  - a list of waves provided by a wave name pattern
      `Xeqt4WList("NameBegin\*NameEnd","Cmd1 ~; Cmd2 ~;")`
  
  - all traces in a graph `Xeqt4TList("Cmd1 §T§")`
  - each integer in an interval `Xeqt4Series(firstNumber,lastNumber,interval,"Cmd1 ~; cmd2 ~;)`
     
 ### *ii) Wrapper function to execute a list of commands in each datafolder inside the current data folder.*
  
  Inside the 
     commands, wildcards, such as, §SUB§ §SUBFULL§, \# stand for the partial or full path to the datafolder,
     the index of the datafolder in the directory etc. 
     
     `XeqtInSubs("Cmd1, Cmd2")`
  
  ### *iii) Functions that add functionality to a graph or automate styling.* 
  - add buttons to scan through all traces, toggling visibility or opacity of traces or pairs of traces. `TraceScanner()`
  - functions to apply custom color schemes with  custom periodicity `ColorByGeo(cyclelength=foo)` `SoftColors()` `GlobalOpacity()`
  - function to display each column of a matrix as a separate trace `PlotTraces()`
  - function to normalize traces by various criteria `NormTraces()`
  - function to automatically distribute multiple vertical axes `DistributeAxes([spacinginpercent])`
  - function to auto-tag traces in a graph `TagTraceWithWaveName()`
  - function to individually color traces according to the values of another wave `ColorTraceByWave([CB])`
  
  ### *iv) Data stratification tools*
  Pick data, either entire waves, or a subset of their entries, or both, based on criteria. This could be numerical entries in the waves notes, or criteria contained in separate waves. Specifically:
  - `WaveSubsetByCriteria(sourceWave, CriteriumWave0, lowCrit0, upCrit0[, CriteriumWave1, lowCrit1, upCrit1, Logic] )` Returns the subset of points in sourceWave for which the corresponding entries in the criteria waves comply with the criteria boundaries
  - `ListDataByNote(Suffix, NoteKey,low4NoteVal, up4NoteVal)` returns waves that have the suffix "Suffix" after a dash, and also have, in their wavenotes the Keyord "NoteKey", followed by ":value" where value lies between low4NoteVal and up4NoteVal
  - `CollectDataByNoteAndCriterium(Suffix[, NoteKey,low4NoteVal, up4NoteVal, CriteriumSuffix, lowCrit, upCrit, CriteriumSuffix1, lowCrit1, upCrit1	] )`combines the above stratification tools in one
  - `QuantilesFromSample()` returns estimates for specified quantiles for a given set of samples.

## Dynamic gain code
This folder contains all functions required to calculate dynamic gain functions from the input and output (current and voltage), to calculate confidence intervals and noise-floor curves for the dynamic gain and to decompose the dynamic gain ([Zhang et al. 2021](https://doi.org/10.1101/2022.02.04.479104)).
Other code provides tools to detect and characterize action potentials (threshold, height, width, depth and time of after-hyperpolarization, peak upstroke and downstroke speed) and action potential time series (coefficient of variation **CV** and local variability **LV** of inter-spike intervals).

## DynamicTimeWarping
This method of waveform-based data alignment is unrelated to dynamic gain calculation. The code in DynamicTimeWarping.ipf is the translation of Matlab code from David Schultz, DAI-Lab, TU Berlin, Germany, 2016, into Igor Pro.
