import pandas as pd
import os
import subprocess
from datetime import datetime, date
import json
import argparse

'''
If it hasn't been done already - create a conda environment with the provided .yaml file and activate the environment.

== Variable dictionary ==

terraWorkspace - input the corresponding Terra workspace where the workflow submissions will be sent, in the format of {NAMESPACE}/{WORKSPACE}, 
this can be found at the top of the page
terraBucket - the bucket address of the corresponding Terra space.
masterSamplesheet - the core sheet that will supply necessary information for the script to fill in and replace, then run the pipeline based on the 
information in the sheet.

== Additional notes ==

For the initial run, select which samples to run via run_pipeline, and use the appropriate sheet (filling in the link_id column with the sampleIDs for multiome.)
Fill in Lane (if applicable), SI_Index column (if applicable, otherwise fill in the indices), project, method, submethod, reference, chemistry (auto if unknown),
flowcell, and path of the BCL directory.

General usage includes 'python3 scRNA_workflow.py [bclconvert | counts | cellbender | cumulus]' (optional) --date  (optional) --submit
If using 'counts', supply --multiome if there are corresponding GEX and ATAC samples, the default option will assume there is no linkage.
If using 'cellbender', supply --post_arc to point the input paths to the corresponding cellranger_arc folders, leave it out if it was not paired data.

Steps should be run sequentially, as the outputs from each step feed into each other, but to just generate outputs and not submit the workflow, omit the --submit flag.

To use a previously generated set of samplesheets and scripts specify the date (--date) to select that set of inputs.

'''

#Core settings that don't need to be changed often
terraWorkspace = "namespace(billing-project)/workspace"
terraBucket = "fc-bucket-here"
masterSamplesheet = "./test.csv"
defaultDate = datetime.now().strftime("%Y_%m_%d")


masterDf = pd.read_csv(masterSamplesheet)
masterDf = masterDf[masterDf['run_pipeline'] == True]
if "link_id" in masterDf.columns:
    masterDf['sampleid'] = masterDf['link_id'] + "_" + masterDf['submethod']

projName = list(masterDf['project'])[0]

mergedIndexJson = {}
jsonsList = [x for x in os.listdir("./ref/") if x.endswith(".json")]
for jsons in jsonsList:
    with open(f"./ref/{jsons}", "r") as f:
        indexs = json.load(f)
        mergedIndexJson.update(indexs)

#Take arguments from script submission and pass into individual steps

def parse_args():
    parserGlobal = argparse.ArgumentParser(add_help=False)
    parserGlobal.add_argument("--date", 
                              default=defaultDate, 
                              help="Enter the date of the folder in YYYY_MM_DD to run the settings from that date, the default dated folder will be created today.")
    
    parserGlobal.add_argument("--submit",
                              required=False,
                              action="store_true",
                              help="Immediately submit run to Terra after creating input and upload inputs")

    parser = argparse.ArgumentParser(description="Running the scRNA workflow one step at a time, each step requires the prior step ",
                                     parents=[parserGlobal])
    subparsers = parser.add_subparsers(dest="step", 
                                       required=True)

    parserBcl = subparsers.add_parser("bclconvert", 
                                      help="Run BCL conversion",
                                      parents=[parserGlobal])
    parserBcl.add_argument("--rev",
                           action="store_true",
                           help = "For use if reverse complement is needed")

    parserCr = subparsers.add_parser("counts",
                                      help="Run Cellranger counts",
                                      parents=[parserGlobal])
    parserCr.add_argument("--modality", 
                          choices=["multiome", "non-multiome", "vdj"], 
                          default="non-multiome", 
                          help="Specify modality for counts")
    parserCr.add_argument("--introns",
                          required = False,
                          action="store_true",
                          help = "If wanting to include introns in cellranger calculations provide this flag.")

    parserCb = subparsers.add_parser("cellbender", 
                                     help="Run Cellbender",
                                     parents=[parserGlobal])
    parserCb.add_argument("--post_arc",
                          action="store_true",
                          help="Add this parameter if the counts results are from Cellranger_arc, this changes the path set to match up with ARC file paths")

    parserCc = subparsers.add_parser("cumulus",
                                     help="Run Cumulus doublet prediction",
                                     parents=[parserGlobal])

    return parser.parse_args()

def preRunCheck():
    if not os.path.exists(f"./scripts/") or not os.path.exists(f"./samplesheets"):
        try:
            os.mkdir(f"./scripts/")
            os.mkdir(f"./samplesheets")
        except:
            print("Error creating samplesheet and scripts folders (try removing them before running the script agains)")
    if not os.path.exists(f"./scripts/{setDate}"):
        os.mkdir(f"./scripts/{setDate}")
    if not os.path.exists(f"./samplesheets/{setDate}"):
        os.mkdir(f"./samplesheets/{setDate}")    

def convertSIToIndex(runDataFrame, indexDict=mergedIndexJson, reverseComplementSeq=False):
    '''
    runDataFrame - Samplesheet that is being processed and read in by the gexAtacSplit() function
    indexDict - the pre-created index dictionary from 10X
    reverseComplementSeq - the deciding factor for whether the other index will be used for RNA/GEX workflows
    '''
    
    def extractIndex(si_index, key=None):
        entry = indexDict.get(si_index)
        if entry is None:
            return None
        if isinstance(entry, dict):   # GEX
            return entry.get(key, None)
        elif isinstance(entry, list): # ATAC
            return entry              # return whole list
        else:
            return None

    for i, row in runDataFrame.iterrows():
        si = row.get("SI_Index")

        # Skip if SI_Index is missing
        if not isinstance(si, str) or si.strip() == "":
            continue

        # --- ATAC ---
        if si.startswith("SI-NA"):
            indices = indexDict.get(si, [])
            if len(indices) != 4:
                raise ValueError(f"Expected 4 indices for {si}, got {indices}")

            runDataFrame.at[i, 'index']  = indices[0]
            runDataFrame.at[i, 'index2'] = indices[1]
            runDataFrame.at[i, 'index3'] = indices[2]
            runDataFrame.at[i, 'index4'] = indices[3]

        # --- GEX ---
        else:
            runDataFrame.at[i, 'index'] = extractIndex(si, "index(i7)")
            if reverseComplementSeq:
                runDataFrame.at[i, 'index2'] = extractIndex(si, "index2_workflow_b(i5)")
            else:
                runDataFrame.at[i, 'index2'] = extractIndex(si, "index2_workflow_a(i5)")

    return runDataFrame

def gexAtacSplit(masterDf):
    '''
    masterDf - main samplesheet as read in as a dataframe
    '''
    dfWithIndexes = convertSIToIndex(masterDf)
    dfWithIndexes.to_csv(masterSamplesheet, index=False)

    createdSheets = {}
    for flowcellId, flowcellDf in dfWithIndexes.groupby("flowcell"):
        
        # --- Decide ATAC vs GEX ---
        def is_atac(row):
            # If SI_Index is present, use it
            if pd.notna(row.get("SI_Index")) and str(row["SI_Index"]).startswith("SI-NA"):
                return True
            # Otherwise: if SI_Index is missing, but index/index2/etc. exist → treat as GEX by default
            return False

        # Split using the helper
        atacDf = flowcellDf[flowcellDf["submethod"].str.upper() == "ATAC"].copy()
        gexDf  = flowcellDf[flowcellDf["submethod"].str.upper() == "RNA"].copy()

        def decideLaneMode(sub_df):
            return (
                "Lane" in sub_df.columns
                and sub_df["Lane"].notna().any()
                and (sub_df["Lane"] != "*").any()
            )

        # --- ATAC ---
        if not atacDf.empty:
            laneModeAtac = decideLaneMode(atacDf)
            if laneModeAtac:
                atacToRun = atacDf[["Lane", "sampleid", "index", "index2", "index3", "index4"]]
                atacToRunL = f"{flowcellId}_ATAC_L_samplesheet.csv"
                atacToRun.to_csv(f"./samplesheets/{setDate}/{atacToRunL}", index=False)
            else:
                atacToRun = atacDf[["sampleid", "index", "index2", "index3", "index4"]]
                atacToRunL = f"{flowcellId}_ATAC_samplesheet.csv"
                atacToRun.to_csv(f"./samplesheets/{setDate}/{atacToRunL}", index=False)
            createdSheets[flowcellId] = atacToRunL
        # --- GEX ---
        if not gexDf.empty:
            laneModeGex = decideLaneMode(gexDf)
            if laneModeGex:
                gexToRun = gexDf[["Lane", "sampleid", "index", "index2"]]
                gexToRunL = f"{flowcellId}_GEX_L_samplesheet.csv"
                gexToRun.to_csv(f"./samplesheets/{setDate}/{flowcellId}_GEX_L_samplesheet.csv", index=False)
            else:
                gexToRun = gexDf[["sampleid", "index", "index2"]]
                getToRunL = f"{flowcellId}_GEX_samplesheet.csv"
                gexToRun.to_csv(f"./samplesheets/{setDate}/{flowcellId}_GEX_samplesheet.csv", index=False)
            createdSheets[flowcellId] = gexToRunL
            
    return createdSheets

def appendSamplesToTemplate(sampleCsv, outputCsv, sep=",", templatePath="./templates/template_samplesheet.csv"):
    '''
    sampleCsv - per sample csv, read in previous output of gexAtacSplit()
    outputCsv - samplesheet names generated as a result of gexAtacSplit()
    sep - comma seperation indicator for pandas
    templatePath - default template path for the samplesheet to be read in
    '''
    
    df = pd.read_csv(sampleCsv)
    rename_map = {
        'lane': 'Lane',
        'sampleid': 'Sample_ID',
    }
    df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns}, inplace=True)

    # Decide order + header string
    if "ATAC" in sampleCsv:
        if 'Lane' in df.columns:
            df = df[['Lane', 'Sample_ID', 'index', 'index2', 'index3', 'index4']]
            header_str = "Lane,Sample_ID,index"
        else:
            df = df[['Sample_ID', 'index', 'index2', 'index3', 'index4']]
            header_str = "Sample_ID,index"

        block1 = df[["Sample_ID", "index"]]
        df = df.drop(columns=["index"])
        block2 = df.melt(
            id_vars=["Sample_ID"],
            value_vars=["index2", "index3", "index4"],
            var_name="index_type",
            value_name="index"
        )[["Sample_ID", "index"]]
        df = pd.concat([block1, block2], ignore_index=True)

    else:
        if 'Lane' in df.columns:
            df = df[['Lane', 'Sample_ID', 'index', 'index2']]
            header_str = "Lane,Sample_ID,index,index2"
        else:
            df = df[['Sample_ID', 'index', 'index2']]
            header_str = "Sample_ID,index,index2"

    # Read template
    with open(templatePath, "r") as f:
        originalLines = f.readlines()
    
    # Custom ATAC block (with spacing)
    atacReadsBlock = [
        "[Reads]\n",
        "Read1Cycles,50\n",
        "Read2Cycles,49\n",
        "Index1Cycles,8\n",
        "Index2Cycles,24\n",
        "\n",
        "[BCLConvert_Settings]\n",
        "CreateFastqForIndexReads,1\n",
        "TrimUMI,0\n",
        "OverrideCycles,Y50;I8;U24;Y49\n"
    ]
    
    with open(outputCsv, "w") as out:
        skipBlock = False
        for line in originalLines:
            stripped = line.strip()

            # Insert ATAC block right after FileFormatVersion
            if "ATAC" in sampleCsv and stripped.startswith("FileFormatVersion"):
                out.write(line)  # keep FileFormatVersion line
                out.write("\n")
                out.writelines(atacReadsBlock)
                continue

            # If we hit [BCLConvert_Settings] in template → skip whole block
            if "ATAC" in sampleCsv and stripped.startswith("[BCLConvert_Settings]"):
                skipBlock = True
                continue

            if skipBlock:
                # End skipping once we reach a new section (starts with "[")
                if stripped.startswith("["):
                    skipBlock = False
                else:
                    continue

            # Stop before [BCLConvert_Data]
            if stripped.startswith("[BCLConvert_Data]"):
                break 

            # Otherwise, copy line
            out.write(line)

        # Now write data section
        out.write("[BCLConvert_Data],\n")
        out.write(header_str + "\n")
        df.to_csv(out, index=False, header=False, sep=sep)

    return outputCsv

#Requires that gexAtacSplit() is run
def configSetupBCL(batchOnly, samplesheets, fc_bucket=terraBucket, runSteps=False):
    """
    batchOnly - DataFrame subset for this batch (must contain 'seq_dir')
    samplesheets - Should be the output of gexAtacSplit in a dict form from the return/result of the previous function
    fc_bucket - workspace bucket
    """
    fcBucket = terraBucket
    bclPaths = set(list(batchOnly['seq_dir']))  # unique flowcell dirs
    configsOut = []

    with open("./templates/bcl_convert_template.json", "r") as f:
        bclJson = json.load(f)

    # Iterate over flowcells defined in samplesheets
    for flowcellId, sheetName in samplesheets.items():
        # Find matching BCL path(s) containing this flowcell ID
        fcMapDf = pd.read_csv(f"./samplesheets/{setDate}/{sheetName}")
        if "Lane" in fcMapDf.columns:
            laneSplit = True
        else:
            laneSplit = False
            
        matching_paths = [p for p in bclPaths if flowcellId in p]
        print(matching_paths)
        if not matching_paths:
            print(f"No matching BCL path found for {flowcellId}, skipping...")
            continue

        for bclPath in matching_paths:
            newConfig = dict(bclJson)  # base config
            newConfig['bclconvert.input_bcl_directory'] = bclPath
            newConfig['bclconvert.no_lane_splitting'] = laneSplit
            newConfig['bclconvert.output_directory'] = f"gs://{fcBucket}/{projName}/fastqs_{projName}/"
            newConfig['bclconvert.sample_sheet'] = f"gs://{fcBucket}/{projName}/samplesheets_{projName}/BCL_Convert_{sheetName}"

            # Write config with flowcell ID in name
            outFile = f"BCLConvert_{flowcellId}.json"
            with open(f"./scripts/{setDate}/{outFile}", "w") as f:
                json.dump(newConfig, f, indent=2)

            configsOut.append(outFile)

            # Optionally push the sheet up to GCP
            if runSteps:
                subprocess.run(
                    f"gcloud storage cp ./samplesheets/{setDate}/BCL_Convert_{sheetName} gs://{fcBucket}/{projName}/samplesheets_{projName}/",
                    shell=True
                )

    return configsOut

def createCountsSheet(countsMode="non-multiome", runSteps=False):
    '''
    countsMode - set by the flag/parameter in the command line submission, changes the input and whether cellranger or cellranger_arc is run
    runSteps - set by --submit, will upload the required files for the cellranger workflow to run on Terra.
    '''
    fcOutputDict = {}
    for fc in set(list(masterDf['flowcell'])):
        with open(f"./scripts/{setDate}/BCLConvert_{fc}.json", "r") as f:
            bclJson = json.load(f)

        fcOutputDict[fc] = bclJson['bclconvert.output_directory'] 

    masterDf['Counts_Input'] = (
        masterDf['flowcell'].map(fcOutputDict)
        + masterDf['seq_dir'].str.rsplit('/', n=1, expand=True)[1]
        + "_fastqs/sample_fastqs/"
        + masterDf['sampleid']
    )

    if countsMode == "non-multiome":
        countsDf = masterDf[masterDf['submethod'] == "rna"]
        countsDf = countsDf[['sampleid', 'Counts_Input', 'chemistry', 'submethod', 'reference']]
        countsDf = countsDf.rename(columns={'sampleid' : 'Sample',
                                            'Counts_Input' : 'Flowcell',
                                            'chemistry' : 'Chemistry',
                                            'submethod' : 'DataType',
                                            'reference' : 'Reference'})
        countsDf.to_csv(f"./samplesheets/{setDate}/cellranger.csv", index=False)
        countsSSLocation = f"gs://{terraBucket}/{projName}/processed/cellranger.csv"
        if runSteps: 
            subprocess.run(f"gcloud storage cp ./samplesheets/{setDate}/cellranger.csv gs://{terraBucket}/{projName}/processed/", shell=True)
            
    elif countsMode == "multiome":
        countsDf = masterDf[['sampleid', 'Counts_Input', 'chemistry', 'submethod', 'reference', 'link_id']]
        countsDf = countsDf.rename(columns={'sampleid' : 'Sample',
                                            'Counts_Input' : 'Flowcell',
                                            'chemistry' : 'Chemistry',
                                            'submethod' : 'DataType',
                                            'reference' : 'Reference',
                                            'link_id' : 'Link'})
        countsDf.to_csv(f"./samplesheets/{setDate}/cellranger_arc.csv", index=False)
        countsSSLocation = f"gs://{terraBucket}/{projName}/processed/cellranger_arc.csv"
        if runSteps:
            subprocess.run(f"gcloud storage cp ./samplesheets/{setDate}/cellranger_arc.csv gs://{terraBucket}/{projName}/processed/", shell=True)
    
    elif countsMode == "vdj":
        vdjMapping = {
            'bcr' : 'vdj_b',
            'tcr' : 'vdj_t'
        }
        masterDf['DataType'] = masterDf['submethod'].apply(vdjMapping)
        countsDf = masterDf[(masterDf['submethod'] == "tcr" or masterDf['submethod'] == "bcr")]
        countsDf = countsDf['sampleid', 'Counts_Input', 'chemistry', 'DataType', 'reference']
        countsDf['Chemistry'] = "fiveprime"
        countsDf = countsDf.rename(columns={'sampleid' : 'Sample',
                                            'Counts_Input' : 'Flowcell',
                                            'chemistry' : 'Chemistry',
                                            'reference' : 'Reference'})
        countsDf.to_csv(f"./samplesheets/{setDate}/cellranger_vdj.csv", index=False)
        countsSSLocation = f"gs://{terraBucket}/{projName}/processed/cellranger_vdj.csv"
        if runSteps:
            subprocess.run(f"gcloud storage cp ./samplesheets/{setDate}/cellranger_vdj.csv gs://{terraBucket}/{projName}/processed/", shell=True)
                           
    return countsSSLocation

def createCumulusSheet(samplesToRun, fcBucket=terraBucket, runSteps=False):
    '''
    samplesToRun - the main samplesheet read in as a dataframe, subsetted to only the RNA/GEX samples due to typical cellbender workflows
    fcBucket - workspace bucket
    runSteps - submitting the job to Terra, controlled by --submit, default function to not submit
    '''
    cbPath = f"gs://{fcBucket}/{projName}/processed/cellbender_v3_{projName}/"
    ccDf = samplesToRun[samplesToRun['submethod'] == 'rna']
    ccDf['Location'] = cbPath + ccDf['sampleid']+ "/" + ccDf['sampleid']+ "_out_filtered.h5"
    ccDf = ccDf[['sampleid', 'Location']]
    ccDf = ccDf.rename(columns={'sampleid' : "Sample"})
    
    for sampleId in ccDf['Sample']:
        tempDf = ccDf[ccDf['Sample'] == sampleId]
        tempDf.to_csv(f"./samplesheets/{setDate}/Cumulus_{sampleId}.csv", index=False)
    
        if runSteps:
            subprocess.run(f"gcloud storage cp ./samplesheets/{setDate}/Cumulus_{sampleId}.csv gs://{fcBucket}/{projName}/processed/cellbender_cumulus_{projName}/{sampleId}/Cumulus_{sampleId}.csv", shell=True)

def configSetupCounts(countsMode, countsSampleSheetPath, intronStatus = "false", fcBucket=terraBucket):
    '''
    countsMode - set by the flag/parameter in the command line submission, changes the input and whether cellranger or cellranger_arc is run
    countsSampleSheetPath - string provided from createCountsSheets() return output, location the samplesheet was uploaded to
    intronStatus - set by the flag/parameter, defaulted to false
    fcBucket - workspace bucket
    '''

    with open("./templates/cellranger_template.json", "r") as f:
        crJson = json.load(f)
    
    if countsMode == "multiome":
        outputLoc = f"gs://{fcBucket}/{projName}/processed/cellranger_arc_{projName}/"
        cellRangerJsonName = "Cellranger_arc"
    elif countsMode == "non-multiome":
        outputLoc = f"gs://{fcBucket}/{projName}/processed/cellranger_{projName}/"
        cellRangerJsonName = "Cellranger"
    elif countsMode == "vdj":
        outputLoc = f"gs://{fcBucket}/{projName}/processed/cellranger_vdj_{projName}/"
        cellRangerJsonName = "Cellranger_vdj"

    newCRJson = dict(crJson)
    newCRJson['cellranger_workflow.input_csv_file'] = countsSampleSheetPath
    newCRJson['cellranger_workflow.output_directory'] = outputLoc
    newCRJson['cellranger_workflow.include_introns'] = intronStatus
        
    with open(f"./scripts/{setDate}/{cellRangerJsonName}.json", "w") as f:
        json.dump(newCRJson, f, indent=2)
    
    crJsonsGenerated = [f"{cellRangerJsonName}.json"]
    
    return crJsonsGenerated

def configSetupCellbender(samplesToRun, postCellrangerArc=False, fcBucket=terraBucket):
    '''
    samplesToRun - the main samplesheet read in as a dataframe, subsetted to only the RNA/GEX samples due to typical cellbender workflows
    postCellrangerArc - a control for directing the path variables to the corresponding "raw_feature_bc_matrix.h5" as an out from standard Cellranger or Cellranger_arc
    fcBucket - workspace bucket
    '''
    with open("./templates/cellbender_template.json") as f:
        cbJson = json.load(f)
    
    samplesToRun = samplesToRun[samplesToRun['submethod'] == 'rna']
    
    cbJsonsGenerated=[]
    for _ , row in samplesToRun.iterrows():
        config = dict(cbJson)
        config["cellbender_remove_background.run_cellbender_remove_background_gpu.sample_name"] = row['sampleid']
        config["cellbender_remove_background.run_cellbender_remove_background_gpu.output_bucket_base_directory"] = f"gs://{fcBucket}/{projName}/processed/cellbender_v3_{projName}/"
        if postCellrangerArc:
            config["cellbender_remove_background.run_cellbender_remove_background_gpu.input_file_unfiltered"] = f"gs://{fcBucket}/{projName}/processed/cellranger_arc_{projName}/{row['link_id']}/raw_feature_bc_matrix.h5"
            outputFile = f"Cellbender_{row['link_id']}.json"
        else:
            config["cellbender_remove_background.run_cellbender_remove_background_gpu.input_file_unfiltered"] = f"gs://{fcBucket}/{projName}/processed/cellranger_{projName}/{row['sampleid']}/raw_feature_bc_matrix.h5"
            outputFile = f"Cellbender_{row['sampleid']}.json"
        
        with open(f"./scripts/{setDate}/{outputFile}", "w") as f:
            json.dump(config, f, indent=2)
        cbJsonsGenerated.append(outputFile)
        
    return cbJsonsGenerated

def configSetupCumulus(samplesToRun, fcBucket=terraBucket):
    '''
    samplesToRun - the main samplesheet read in as a dataframe, subsetted to only the RNA/GEX samples due to typical cellbender workflows
    postCellrangerArc - a control for directing the path variables to the corresponding "raw_feature_bc_matrix.h5" as an out from standard Cellranger or Cellranger_arc
    fcBucket - workspace bucket
    '''

    with open("./templates/cumulus_template.json") as f:
        ccJson = json.load(f)
    
    samplesToRun = samplesToRun[samplesToRun['submethod'] == 'rna']

    ccJsonsGenerated = []

    for _, row in samplesToRun.iterrows():
        config = dict(ccJson)
        config['cumulus.input_file'] = f"gs://{fcBucket}/{projName}/processed/cellbender_cumulus_{projName}/{row['sampleid']}/Cumulus_{row['sampleid']}.csv"
        config['cumulus.output_directory'] = f"gs://{fcBucket}/{projName}/processed/cellbender_cumulus_{projName}/"
        config['cumulus.output_name'] = f"{row['sampleid']}"
        outputFile = f"Cumulus_{row['sampleid']}.json"
        with open(f"./scripts/{setDate}/{outputFile}", "w") as f:
            json.dump(config, f, indent=2)
        ccJsonsGenerated.append(outputFile)

    return ccJsonsGenerated

def scriptSteps(step, jsonGenerated, runStep=False):
    '''
    step - steps fed in from argparse, of three steps (see additional notes on github)
    jsonGenerated - outputs from the ([configSetupBCL, configSetupCounts, configSetupCellbender]) functions, fed directly due to structuring for the alto terra submission script.
    runStep - submitting the job to Terra, controlled by --submit, default function to not submit
    '''
    stepDict = {"bclconvert" : "kco/bcl_convert",
                "counts" : "lilab:cumulus:Cellranger:3.1.1",
                "cellbender" : "cellbender/remove-background/13",
                "cumulus" : "cumulus/cumulus"}
    
    scriptBaseText = f'''#!/bin/bash
    
alto terra run \
    -w {terraWorkspace} \
    -m {stepDict[step]} \
    -i JsonValue
'''
    
    scriptsGenerated = []
    for scriptJsonValue in jsonGenerated:
        scriptBaseTextEdit = scriptBaseText.replace("JsonValue", f"./scripts/{setDate}/{scriptJsonValue}")

        completePathScripts = f'./scripts/{setDate}/alto_terra_{scriptJsonValue.replace(".json",".sh")}'
        scriptsGenerated.append(completePathScripts)
        
        with open(completePathScripts, 'w') as f:
            f.write(scriptBaseTextEdit)
    
    if runStep:
        for scripts in scriptsGenerated:
            subprocess.run(f"bash {scripts}", shell=True)
    
    return scriptsGenerated

if __name__ == "__main__":
    args = parse_args()

    if args.date != defaultDate:
        setDate = args.date
    else:
        setDate = defaultDate

    print(f"Date for run: {args.date}")

    preRunCheck()
    print(f"Step: {args.step}")

    if args.step == "bclconvert":
        convertSIToIndex(masterDf, mergedIndexJson, args.rev)
        createdSheets = gexAtacSplit(masterDf)
        for sheets in createdSheets.values():
            appendSamplesToTemplate(f"./samplesheets/{setDate}/{sheets}",
                                    f"./samplesheets/{setDate}/BCL_Convert_{sheets}",
                                    templatePath="./templates/template_samplesheet.csv")
        jsonGen = configSetupBCL(masterDf, createdSheets, runSteps=args.submit)
    
    elif args.step == "counts":
        print(f"Counts modality: {args.modality}")
        countsSamplesLocation= createCountsSheet(args.modality, args.submit)
        jsonGen = configSetupCounts(countsMode=args.modality,
                                    intronStatus = args.introns,
                                    countsSampleSheetPath=countsSamplesLocation)

    elif args.step == "cellbender":
        jsonGen = configSetupCellbender(masterDf, args.post_arc)
    
    elif args.step == "cumulus":
        createCumulusSheet(masterDf, runSteps = args.submit)
        jsonGen = configSetupCumulus(masterDf)
    
    #Catch function to have jsonGen be a list when submitted, else it will iterate through a string
    if isinstance(jsonGen, str):
        jsonGen = [jsonGen]
    
    scriptSteps(args.step, jsonGen, args.submit)