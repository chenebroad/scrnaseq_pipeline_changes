This is rewrite for the original scRNA pipeline constructed by Daniel Chafamo, with a less reliance on provisioning virtual machines and a goal of removing reliance on dsub due to changes in Google Batch / Life Sciences API which have made maintaining the original pipeline a little more challenging.

The main script allows for the creation and upload of inputs to the corresponding Terra workspace for single cell processing. 

In principal, the pipeline/script runs similarly where in the workflows are submitted to Terra for parallelization and actual computation, and the script functions as a liasion for creating the inputs and stringing them together for stepwise submission.
A .yaml file has been provided with corresponding packages required to run the script. 

## Installation

Outside of installing the yaml, one would also need to create an environment:

`conda create -n alto-terra python=3.10`

Additionally, making sure that permissions are set up accordingly for the Google Cloud SDK:

`gcloud auth login`

For future considerations, it may be important to include more workflow specific parameters to slot into the WDLs/workflows maintained by the Cumulus and Cellbender teams.

Step-wise there are three main steps that run off of the csv file provided - depending on the incoming data it is imperative to use the correct samplesheet template:
The sample tracking file, in csv format, is a useful way to track the important information for each sample, and is needed to run this script. Each sample requires the following text fields.

## GENERAL TRACKING:

- date (OPTIONAL): The date your samples are processed in yyyy_mm_dd format.
- sampleid (REQUIRED): This is the sample id.
- tissue (OPTIONAL): The tissue of origin.
- replicate (OPTIONAL): No longer used but potentially useful for tracking replicates
- condition (OPTIONAL): No longer used but potentially useful for tracking special collection conditions

Detailed below are columns for each step.

## ALL:

- run_pipeline (REQUIRED): Boolean (True or False) that determines what samples are processed. Set this to True for all samples you want to processs. All other samples must be set to False. How this works in operation is that as you add your new samples, set them to run_pipeline = True and set the previously run samples to run_pipeline = False. Remember that all samples that are processed together must come from the same flow cell. The code is written to only process one flow cell!
project (REQUIRED): The name of the project you'd like to see attached to your directories

## BCL_CONVERT:

- link_id (REQUIRED): Base sample name to be used for linking multiome samples together for post processing, GEX and ATAC samples that are from the same subject should have the sample link_id.
- sample_id (REQUIRED): General name assignment for the samples given to generated downstream files, in multiome samples this will function slightly differently in order to fit accomodations for the cellranger_arc workflow, will serve as a way to uniquely identify linked/paired samples.
- method (OPTIONAL): [rna, atac, or multiome] currently does not control any core functions but will help with delineating any files.
- submethod: [rna, atac] now pridominantly used for determining how to run the files through BCLConvert, where in all indices from ATAC data will be captured and structured accordingly to BCL_CONVERT samplesheet requirements. rna corresponds to GEX data, whilst atac will be for ATAC data.
- index[2,3,4] (CONTINGENT): The corresponding indices written in the nucleotide form.
- flowcell (REQUIRED): The flowcell id from your sequencing run, has slight functionings with path generation and seperation, it is more important to label samples that go together than for actual running processes, serves as a grouping measure. Generally speaking, if a flowcell is given (250825_VH00997_445_222FMY3NX), it will should be the last portion after the last underscore (e.g. 222FMY3NX).
- seq_dir (REQUIRED):  The directory of your sequencing results in GCP, downstream folder names will pull from this varible.

## COUNTS:

See above: sample_id

- reference (REQUIRED): The genome reference to use when Cell Ranger count is creating the counts matrices. Please choose from one of references listed in Cumulus read the docs.
- chemistry (CONTINGENT): The sequencing chemistry used.
- introns (OPTIONAL): Include introns in the cellranger workflow; default set to false unless specified.

## CELLBENDER:

See_above: sample_id

Legacy columns with potential use (pre-Cellbender v0.3.0 required manual assignment of parameters, now it is automaticall calculated.) If need for manual parameter selection is required, the functionality will be added to the current script.

- min_umis (OPTIONAL): the min number of UMIs you'd like to use for filtering when you run cumulus pegasus. 
- min_genes (OPTIONAL): the min number of genes you'd like to use for filtering when you run cumulus pegasus.
- percent_mito (OPTIONAL): the max percentage of expression coming from mito genes that you'd like to set for filtering when you run cumulus pegasus.
- cellbender_expected_cells (OPTIONAL): input of expected cells for Cellbender's cell calling
- cellbender_total_droplets_included (OPTIONAL): input of empty droplets for Cellbender's cell calling
