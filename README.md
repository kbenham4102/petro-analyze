# petro-analyze
Tools for Formula Identification Petroleomics Data Analysis
# Getting Started
Requirements:
<br> - Anaconda 3 package - everything you need should be here.

## Running the Command Line Petroleomics Plot Generator

<br> There are a few arguments to consider for these plots.
<br> To see all arguments in syntax run ```$ python compound_id.py -h``` for the help menu.
<br> Here is the syntax to generate all plots for the example datafiles given:
```
python compound_id.py -data YourDataFile.csv --w 1 --err_plots 1 --dbe_plots 1 --class_plots 1 --keep_class HC,N,O,O2 --comp_data ComparisonDataFile.csv --comp_keep_class HC,N,O,O2 


```
## Note On Class Distribution Plots
The class distribution plots are generated based on a standard list ```class_list``` in the code. You can use the interactive notebook to figure out which classes are most abundant and edit the list manually. This could be an added feature in the fututure to automatically generate plots for standard lists and most abundant classes.

## For Development...
Also included is the ipython notebook ```compound_id.ipynb``` where you can develop new plots with interactively in Jupyter. Examples used in the current program can be found here.

## TODO
<br> Van krevelen plots
<br> Kendrick mass plots
