# OESpikeUtilities

This plugin packages provides several plugins used for online spike detection and display. The plugins include

## Dynamic Spike Detector
This plugin was modified from the [original dynamic threshold plugin](https://github.com/camigord/DynamicSpikeDetectorPlugin). We fixed some of its bugs and added in various quality of life improvements. In particular, it allows setting the threshold for all channel at once. It also allow specifying the detection sign. The spike detection only works when the `Enable` toggle is enabled. It currently supports two methods of spike detection. The `Simple` one just detect spike based on signal ampitude. The `Median` method is based on median of the absolute signal. Please note the `Simple` one may actually work better on whitened signals. The `Median` method calculate the median threshold based on a very small segment of signal in the input buffer, which may not be optimum. We plan to change this behaviour in near future.

<img width="333" alt="image" src="https://user-images.githubusercontent.com/3406709/146046735-ac352dac-4223-4c53-abbe-a65ad93d523f.png">

## LFP Spike Viewer

This plugin is modified from the original LFP viewer but add an overlay on top of the neural signals to show the spike location. It will be useful for monitoring the quality of the spike detection and help setting the detect threshold.
<img width="914" alt="image" src="https://user-images.githubusercontent.com/3406709/146047977-8ea3aab3-73c8-47cf-be67-8cdf9aee85b7.png">

## Installation
Please follow the Open Ephys GUI page to build this plugin: [link](https://open-ephys.atlassian.net/wiki/spaces/OEW/pages/1301643269/Creating+Build+files)
