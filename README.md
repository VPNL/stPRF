# Spatiotemporal pRF 

Contact: Insub Kim (insubkim@stanford.edu)

## Dependencies
Scripts were developed using MATLAB 2020b (Mac OS/Linux).

Toolboxes: [SpatiotemporalpRFs](https://github.com/VPNL/spatiotemporalPRFs), [PRFmodels](https://github.com/vistalab/PRFmodel), [Vistasoft](https://github.com/KimInsub/vistasoft.git), [BADS](https://github.com/acerbilab/bads.git)

You can download the all the toolboxes by running:
```Matlab
% MATLAB: download_toolboxes.m
download_toolboxes('toolbox_urls.txt', 'code/toolbox')
```
It is recommended to run it on a machine with a **GPU** to solve the spatiotemporal pRF estimates


## Experiment
[![experiment](https://img.youtube.com/vi/M9LGkSw4n8Q/0.jpg)](https://youtu.be/M9LGkSw4n8Q)  
Click to view it on YouTube

## Implemented pRF models
<img src="doc/models.png " width="700">

## Demo

```Matlab
%% spatiotemporal pRF Demo (run.m)
% - download toolboxes
% - set variables and create JSON files
% - create synthetic timecourses
% - solve spatiotemporal PRF models
% - plot results

run.m

```

### (1) Synthetic timecourse generation
The software takes stimulus information and a JSON file ([example JSONfile](https://github.com/VPNL/stPRF/blob/main/template/template.json)) as input and generates synthetic timecourses with noise for the three different pRF models: `spatial`, `DN-ST`, and `CST`.

**JSON file:** To compute spatiotemporal pRF models, the JSON file's temporal class must be specified. 
Three different types of pRF models can be used: spatial pRF model, CST model, and DN-ST model. The spatial pRF model is similar to the one-Gaussian pRF model but uses a higher temporal sampling resolution. CST and DN-ST models are spatiotemporal pRF models that have explicit temporal impulse response functions. 

For more detailed information about how to edit the JSON file ([Link](https://github.com/vistalab/PRFmodel/wiki/prf-Synthesize:-how-to-edit-json-file))

### (2) Solve pRF models
Solve the parameters for each model(`spatial`, `DN-ST`, and `CST`). Synthetic timecourses generated for each model are solved by the same model.

### (3) Check performance
Compare the ground truth and predicted timecourses for each model, and plot the results.


## Paper

* Code to regenerate figures in the paper [Link](https://github.com/VPNL/stPRF/tree/main/code/figures)

## References
Dumoulin, S. O., & Wandell, B. A. (2008). Population receptive field estimates in human visual cortex. Neuroimage, 39(2), 647-660.

Lerma-Usabiaga, G., Benson, N., Winawer, J., & Wandell, B. A. (2020). A validation framework for neuroimaging software: The case of population receptive fields. PLoS computational biology, 16(6), e1007924.

Zhou, J., Benson, N. C., Kay, K., & Winawer, J. (2019). Predicting neuronal dynamics with a delayed gain control model. PLoS computational biology, 15(11), e1007484.

Stigliani, A., Jeska, B., & Grill-Spector, K. (2019). Differential sustained and transient temporal processing across visual streams. PLoS computational biology, 15(5), e1007011.
