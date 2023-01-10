# Oxygen_Inverse_Model

The following code implements an inverse model that estimates d18O of body water from d18O measurement series obtained from sequentially sampled tooth enamel, which was originally published in Passey, Benjamin H., et al. "Inverse methods for estimating primary input signals from time-averaged isotope 
profiles." Geochimica et Cosmochimica Acta 69.16 (2005): 4101-4116.https://doi.org/10.1016/j.gca.2004.12.002

This original publication implemented the inverse model using MATLAB, which is relatively expensive proprietary software. The code in this repository represents an adaptation of the model as is into R. It draws heavily on the R translation of the original MATLAB code published in Fraser, Danielle, et al. "Pronghorn (Antilocapra americana) enamel phosphate δ18O values reflect climate seasonality: Implications for paleoclimate reconstruction." Ecology and Evolution 11.23 (2021): 17005-17021. https://doi.org/10.1002/ece3.8337, but adds some additional features such as generating automatic html reports of model parameters and outcomes and customized plots of model results. 

Variable and parameter naming follows that of the original Passey et al. (2005) code. The exact modelling procedure and what each parameter means is described in detail in the Supplementary Information
for that paper,mostly in Electronic Annex 4. While the code in this repository does contain comments and explanations it is highly advisable to study the model publication and Supplementary Information in detail before using these scripts. 


## Prerequisites

In order to be able to run the code, you will need to have an up-to-date version of [R](https://www.r-project.org/) and [RStudio](https://rstudio.com/) or [pandoc](https://pandoc.org/) installed on your computer and a few CRAN packages (see below). 

The code was last run successfully using R version 4.2.0 (2022-04-22 ucrt) on a 64 bit Windows 10 operating system. 

### Install packages

```
install.packages(c("ggplot2", "dplyr", "matrixcalc", "rmarkdown"))
```

Last successfully used package versions were:

ggplot2_3.3.6    
dplyr_1.0.9      
matrixcalc_1.0-5
rmarkdown_2.16

## Files and code

This code is designed to run in an RStudio project environment and file paths are given relative to the project directory. 

`run_inverse_model.r`

 - R script to execute the inverse modelling procedure. Please refer to the original publications by Passey et al. (2005) and Fraser et al. (2021) for detailed explanations of model equations and parameters. Input parameters need to be changed as needed and the script is run with adjusted parameters once for each specimen of interest. The parameters that were used in this study for each specimen are documented in the model report files contained in `run_inverse_model/output`. 

`make_model_report.Rmd` 

- RMarkdown template called by `run_inverse_model.R`
to print model report html files

`input/tooth1_inverse_model_input.csv`

- example d18Oenamel data to use as input for the inverse model

`functions/Emeas.R`

- R script containing the `Emeasfun` function called by the `run_inverse_model.r` script. This is equivalent to the Emeas portion of the codes provided by Passey et al. (2005)

`functions/mSolv1.R`

- R script containing the `PasseyInverse` function called by the `run_inverse_model.r` script. This is equivalent to the mSolv portion of the codes provided by Passey et al. (2005) and has been adapted from the R translation of the the original Passey Matlab codes that has been published by Fraser et al. (2021). The model portions of this code are identical to those provided by Fraser, but have been updated with increased commenting and modified to export additional model outcomes and enable execution and output processing in external scripts in a more automated way. 

`output/tooth1_inverse_model_output.csv`

- results of the inverse model, including results for all individual trials, best model solution and confidence intervals

`output/tooth1_inverse_model_plot.png`

- plot of the model outcome and reference vector compared to the original d18Oenamel input

`output/tooth1_inverse_model_report.html`

- automated report generated during execution of `run_inverse_model.r`. Documents input parameters and select model diagnostics. 

## Running the inverse model

### Step 1 - make input data

Input data should be formatted the same as the example input provided here. Some columns are not necessary for the model itself, but are useful for contextual information. These are denoted as optional. 

tooth_ID - the ID of the tooth specimen (e.g. find number), used for naming output files

sample_ID - unique sample identifier (optional)

layer - layer information (optional, but code lines referring to this should be removed if this is not provided)

taxon - species (optional)

dist_ERJ - distance from root-enamel-junction in mm

sample_length - the distance between this sample and the previous sample in mm. Computed from dist_ERJ. 

mean_d18O - the oxygen isotope delta value (can be VSMOW or VPDB)

SD_d18O - measurement uncertainty, here SD of replicate measurements

sample_depth - drilling depth either in mm or as fraction of total

enamel depth (the example uses measurements in mm)

### Step 2 - generate Emeas estimates


Use the portion 'running the Emeas function' of the run_model.R script to generate Edist estimates

Start by putting the user input parameters

1) numtrials = 100 # leaving this at 100 is a good idea, more is slow and less is not enough trials
2) set r1, r2, r3 and la paramters correctly depending on input data and taxon specific tooth mineralization (parameters correspond to descriptions in Passey et al.)
3) run the Emeas function (Emeasfun). The output of Emeas values from Emeasfun is stored in an object called 'Edist'. 
4) plot histogram of Edist (analogous to Passey et al. 2005 Electronic Annex 4, Figure EA-4-1) and extract mean of Emeas from Edist. This will be used later to adjust the damping factor df

### Step 3 - adjust damping factor df

Use the portion 'mSolv1_1 script' portion of the run_model.R script. 

In this portion of the script the mSolv function (analagous to Passey mSolv1_1 script) is used to run some initial trials to adjust the damping factor (df).  This portion of the script generates estimates of the prediction Error (Epred), which should be matched to have a similar value to Emeas. The values of Epred are stored in an object called 'DPE'. We therefore run the mSolv function multiple times with a lower number of iterations (~50 to 100) and change df until Emeas (stored in Edist) and Epred (stored in DPE) are close in value. 

1) Adjust mSolv parameters: set number of trials to 50 (nsolxns = 50). Set the reference vector input parameters to be close to the mean of the data with a relatively flat reference vector. Start with a df of 0.01. 

2) run mSolv with 50 solutions
3) A plot of first 10 solutions is used to get a first idea of the shape of the model output. 
4) Check mean of Edist (Emeas values) vs mean of DPE (Epred values). They should be relatively close in value. Also plot here the histogram of DPE, to check how representative the mean is of the distribution. If the DPE values are very scattered it may be necessary to run more trials (e.g. 100), to give more repeatable values for DPE. 
5) Adjust df to match Edist and DPE (higher df means higher DPE). Repeat until Edist and DPE are matched to ~ first decimal. Then proceed to Step 4. 

### Step 4 - check sensitivity to the reference vector


Check the sensitivity to reference vector oscillation (determined by maxratio and minratio. The further they are apart the more oscillation is in the reference vector). It makes sense to check a relatively flat reference vector with ca 1 permil oscillation and then a reference vector with oscillation that is similar to the data. 

1) adjust maxratio and minratio
2) re-run mSolv with the different reference vectors
3) compare the shape of solutions (in blue) in the plot mSolv output between the model runs
4) repeat 1-2 with less/more oscillating ref vector and check for major changes in amplitude. Small changes in solution shape outside the seasonal peaks and troughs can be discounted. 
5) The reference vector should not have major influence on the amplitude of the solutions

### Step 5 - run the full model with adjusted df

1) set nsolxns to at least 200
2) reload all input parameters and empty objects
3) run mSolv function to run model
4) run the 'generate report of model parameters section of the script. This generate an html report of model input parameters and some model diagnostics to document the parameters that were chosen for each d18O time series
5) run the 'extract output with mean and CI and write to csv' portion of the script. This will write model output into a csv file in a form that is usable for further processing (e.g. plotting, extracting values etc)

### Step 6 - plot model results

Run the last portion of the script to plot the result of the inverse model. The plot will display the original stable isotope data, the most likely model solution and model confidence interval as well as the reference vector. 

The result is exported as tooth1_passey_model_plot.png
