# COASTER: Confidence Assessment Tool for eDNA and qPCR Results

(TODO) ADD LICENCE

(TODO) DISCLAIMER ON INTERPRETATION - IMPORTANT

## Overview
The primary purpose of COASTER is to aid in the interpretation of quantitative PCR (qPCR) results from environmental DNA (eDNA) samples. A standardised PDF report will be generated to accompany your qPCR results when sent to other scientists and decision-makers. 

Alongside your qPCR data, COASTER requires two things:


**1. Assay Validation:** COASTER assigns a confidence level (Low|Medium|High) to a qPCR assay based on a simple checklist of the assay validation steps that have been so far undertaken. Broadly, the confidence levels can be viewed as follows:

- Low: interpret both positive and negative results with caution. There can be false positives (a positive signal may have arisen from a non-target species). There can be false negatives (an absence of detection does not indicate that the species is absent from the sampling point).
- Medium: positive samples can interpreted as "Species DNA present in sample". There can be false negatives. 
- High: positive samples can be interpreted as "Species DNA present in sample". There is low risk of false negatives.


**2. qPCR Experimental Settings:** You are then required to fill in simple information on the experimental settings used to run your qPCR plates. 

The script is implemented as a Shiny app which has a publicly available dashboard [here](https://vidasolutions.shinyapps.io/TOOL_dashboard/)


_Note: COASTER remains a work in progress and we will fix new bugs if identified by users of the tool._

## Usage

**Step 1:**
_**Important:**_ The tool requires a very specific data format to function, so please download a modifiable template from the home page. The template consists of a csv file with headers:

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
**Step 2:** Fill out the Assay Validation table by carefully considering and answering each of the questions (1 = Yes, 0 = No).

**Step 3:** Fill out the qPCR Experimental Settings you used to generate your qPCR results.

**Step 4:** If you have correctly downloaded the data template and filled this out locally with the required information (see Appendix 1 for an example), upload this in the 'Data Upload Section' (Figure 2). A table and standard curve plots (for each plate) will then be generated. These can be checked to ensure your data looks as expected.

**Step 5:** Once the upload is complete the user should enter a project title. Additionally, the user may choose to select for the report to include an appendix with interpretation of the results at the qPCR replicate level, which will generate a table with the results inputted to the template file plus an additional column with an interpretation of each replicate result.

**Step 6:** Select generate report. The generated report will give the date and time the report was generated, followed by four sections: summary of results, settings, results and interpretations, and appendix. An example of a generated report can be seen in Appendix 2, but a brief description of each section is as follows:

1.	Summary of results includes the assay validation level, if contamination was observed in any controls, if positive controls amplified as expected, and if standards were included in the uploaded data the associated R2 and PCR efficiency values.
2.	Settings includes a table with the parameters listed as set by the user.
3.	Results and interpretations includes a table summarising the results at 1) the Sampling_Point level 2) the DNA sample level and 3) [optional] the qPCR replicate level 
4.	Samples will either be classified as positive, negative, inconclusive, or tentative. A tentative is the result of weak signal (Cq value below the limit of detection).

