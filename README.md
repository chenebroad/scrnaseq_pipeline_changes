This is rewrite for the original scRNA pipeline constructed by Daniel Chafamo, with a less reliance on provisioning virtual machines and a goal of removing reliance on dsub due to changes in Google Batch / Life Sciences API which have made maintaining the original pipeline a little more challenging.

The main script allows for the creation and upload of inputs to the corresponding Terra workspace for single cell processing. 

In principal, the pipeline/script runs similarly where in the workflows are submitted to Terra for parallelization and actual computation, and the script functions as a liasion for creating the inputs and stringing them together for stepwise submission.
A .yaml file has been provided with corresponding packages required to run the script. 

For future considerations, it may be important to include more workflow specific parameters to slot into the WDLs/workflows maintained by the Cumulus and Cellbender Teams.
