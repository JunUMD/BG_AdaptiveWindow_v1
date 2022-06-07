# BG_AdaptiveWindow_v1
Improving ATMS Remapping Accuracy Using Adaptive Window and Noise-tuning Method in Backus-Gilbert Inversion

This software package is developed to remap ATMS channel 1 observations from 5.2° FOV to a consistent AMSU-A 3.3° FOV. To improve the remapping accuracy, an adaptive window method is applied to provide sufficient information for the reconstruction at each scan position. In addition, a new noise tuning method is proposed to eliminate the scan-angle-dependent features in the noise caused by the sensor’s cross-track scanning manner. Results from simulations and NOAA ATMS data show that compared to the fixed window, the new method can significantly reduce the bias stemming from the resolution difference. The issue of the deterioration of the resolution enhancement capability near the scan edge in the fixed window method has been largely ameliorated.

The code can be easily modified and applied to other spaceborne microwave sensors. Please cite the manuscript “Improving ATMS Remapping Accuracy Using Adaptive Window and Noise-tuning Method in Backus-Gilbert Inversion” published on IEEE TGRS when you use this software package in your own research.
