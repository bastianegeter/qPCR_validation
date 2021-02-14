# eDNA qPCR assay and project confidence assessment tool 

(TODO) ADD LICENCE

(TODO) DISCLAIMER ON INTERPRETATION - IMPORTANT

## Overview
The primary purpose of this tool is to aid in the interpretation of quantitative PCR (qPCR) results from environmental DNA (eDNA) samples. The tool is also used to assign a confidence level (Low, Medium, High) to a qPCR assay based on a simple checklist of the assay validation steps that have been so far undertaken. Broadly, the confidence levels can be viewed as follows:

Low: interpret both positive and negative results with caution. There can be false positives (a positive signal may have arisen from a non-target species). There can be false negatives (an absence of detection does not indicate that the species is absent from the sampling point).
Medium: positive samples can interpreted as "Species DNA present in sample". There can be false negatives. 
High: positive samples can be interpreted as "Species DNA present in sample". There is low risk of false negatives.

The script is implemented as a Shiny app which has a publicly available dashboard [here](https://vidasolutions.shinyapps.io/TOOL_dashboard/)

## Usage
Step 1: A modifiable template can be downloaded from the dashboard. The template consists of a csv file with headers. 

```
Plate	Well	Sample_Type	DNA_Sample	Replicate	Target_Cq	IPC_Cq	Std_Conc	Sampling_Point	Extraction_Batch	Volume_Water_Processed

```

## Explanation of headers
TODO - add specs

```
Plate	The plate ID that the DNA samples are being run on
Well	Location of the sample on the plate
Sample_Type	Type of DNA sample. There are six sample types allowed: external negative control (extnc), PCR Negative control (pcrnc), field negative control (fieldnc), positive control (pc), standard (std), and unknown (unkn)
DNA_Sample	Sample ID number or name
Replicate	Replicate number for a sample
Target_Cq	Cq value obtained for a sample
IPC_Cq	Cq value obtained for the IPC
Std_Conc	Concentration of the standard in ng µl-1
Sampling_Point	Site name or number
Extraction_Batch	Extraction batch number
Volume_Water_Processed	Volume of water filtered for a sample in ml


```

Step 2: After the template has been downloaded and filled out locally with the required information (see Appendix 1 for an example), it can be uploaded in the appropriate section of the homepage (Figure 2).
Step 3: Once the upload is complete the user should enter a project title. Additionally, the user may choose to select for the report to include an appendix with interpretation of the results at the qPCR replicate level, which will generate a table with the results inputted to the template file plus an additional column with an interpretation of each replicate result.
Step 4: The user must now set the appropriate parameters in the ‘Set parameters’ section of the tool homepage (Figure 2). All parameters must have a value entered in the setting column.
Step 5: Select generate report. The generated report will give the date and time the report was generated, followed by four sections: summary of results, settings, results and interpretations, and appendix. An example of a generated report can be seen in Appendix 2, but a brief description of each section is as follows:
1.	Summary of results includes the assay validation level, if contamination was observed in any controls, if positive controls amplified as expected, and if standards were included in the uploaded data the associated R2 and PCR efficiency values.
2.	Settings includes a table with the parameters listed as set by the user.
3.	Results and interpretations includes a table summarising the results at 1) the Sampling_Point level 2) the DNA sample level and 3) [optional] the qPCR replicate level 
4.	Samples will either be classified as positive, negative, inconclusive, or tentative. A tentative is the result of weak signal (Cq value below the limit of detection).

