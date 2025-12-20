<!DOCTYPE html>
<html>
<body>

The Dynamic Gain Code requires the *ANtools_extended.ipf*  but no external functions from outside the repository. Should you wish to import data in the proprietary file formats of *Clampex* or *PatchMaster* into Igor Pro, I suggest you use the additional files created by htasche, which are available at [bpc_ReadAxg](https://www.wavemetrics.com/project/bpc_ReadAxg) and [bpc_ReadHeka](https://www.wavemetrics.com/project/bpc_ReadHeka).

### Requirements

Raw data have to be represented as pairs of voltage and current waves which have to be named with a common trial name and the suffixes "_V" and "_I". For automation, all trial names should start with the same prefix. The units of amplitude and time should be Volt and Ampere, and Seconds.  

All parameters, such as voltage thresholds or minimal time intervals for spike detection are expected to be in SI units (Volts, Seconds, Ampere etc.). The outcome of analyses, spike train features (e.g inter-spike intervals), action potential shape features (e.g. width, height, after-hyperpolarization magnitude), and Dynamic Gain Magnitude and Phase are also given in units that are combinations of SI units without prefactors (e.g. Hertz / Ampere for the Dynamic Gain Magnitude).

### Overview of the implemented calculations
The calculations are explained in the publications, in particular in [Merino et al. 2021 PNAS](https://doi.org/10.1073/pnas.2114549118) and [Zhang et al. 2024 J Neurosci](https://doi.org/10.1523/JNEUROSCI.0799-23.2023). The following short summary gives an overview. More details are provided in the footnotes. A comparison between two different methods are provided in the supporting material in [Revah et al. 2025 bioRxiv](https://www.biorxiv.org/content/10.1101/2024.02.15.580451v5). There, we show that the linear approach captures the response very well. We also show that the Dynamic Gain measurements based on the noise input alone (spike-triggered average current) or based on an additional sinusoidal input (vector strength) fully agree with each other. 


![Main](Flow-chart-Main-broad.png)
![Footnotes](Footnotes.png)

