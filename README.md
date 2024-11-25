Raw data files of the form AXX_pitchX_TimeStamp_TimeStamp.mat should be located in the ./DataTransients/ folder (<root>).

The analysis file Analysis_Transients.m will go through those files and save the results in the folder <root>/processed/ by appending "_PROCESSED" in the filename before the .mat extension.

The analysis file Analysis_PlotProcessedData.m retrieves those _PROCESSED.mat files plus all geometry and similar files in the filesystem tree to plot relevant quantities in scatter plots or boxplots. The first part of the analysis file is just loading and converting units (ending around line 260). The second part can be run cell by cell to produce plots. Each cell's first line contains a short description of what is plotted. Towards the end the first section there is a commented section that explains what variable is located in each dimension of the data matrix MX(..., ..., ...). That matrix contains (100,...,...) nans by default that are later populated with actual data (to avoid defaulting with other quantities that would show in plots); this is done because each data set produces a different amount of valid data, between 40 and 60 transients per dataset. The second dimension of MX(..., 22, ...) corresponds to the 22 dataset files. The third dimension accesses the variable to be plotted.


Required files in the tree include:
In <root>
* Acquisition_Transients.m
* Analysis_PlotProcessedData.m
* Analysis_Transients.m
* DecomposeFourier_func.m
* DetectVelocityPeaks_func.m
* Lp_data.txt
* Lp_m.mat
* myenvelopedetection.m
* Vf_data.txt
* Vf_m3.mat

in <root>/processed 
* Qfac_anal.m (produces the modal maximal quality factors as per Blanc2010 paper)

etc.



## Some description of dependencies
`plot_getScatteredData.m` is used to retrieve all XData and YData from a scatterplot with nans and multiples groups due to MX(transNum, datasetFile, variableNum)

(the fundamental search algorithm 'yin' should be in your MATLAB search path)

## Quick summary ##
`Analysis_Transients.m` will process the raw data and create `_PROCESSED.mat` files.
`Analysis_PlotProcessedData.m` takes the `_PROCESSED.mat` data and aides visualization. It uses `Analysis_PlotProcessedData_Loads.m` for loads and population of variables.
