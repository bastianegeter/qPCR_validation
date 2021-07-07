# COASTER: Confidence Assessment Tool for eDNA and qPCR Results

(TODO) ADD LICENCE

## Overview
Environmental DNA (eDNA) is increasingly being used to survey for local biodiversity. Despite the opportunities that eDNA provides, results can be difficult to interpret and compare between datasets that have used different **_assays_** and **_experimental methods_**. To overcome this, COASTER was developed for scientists and decision-makers to aid in the interpretation of qPCR results from eDNA surveys. The tool assigns a confidence level (Low|Medium|High) to a qPCR **_assay_** based on a simple checklist and automatically interprets the confidence of qPCR results based on the **_experimental methods_** used.

To get started with COASTER, you will need to provide: 

**1. qPCR results:** The qPCR results you would like to be assessed (please see below regarding the _data template_.

**2. Assay Validation:** COASTER assigns a confidence level (Low|Medium|High) to the qPCR assay used based on a simple checklist of the assay validation steps taken so far. Broadly, the confidence levels can be viewed as follows:
- _Low:_ interpret both positive and negative results with caution. There can be false positives (a positive signal may have arisen from a non-target species). There can be false negatives (an absence of detection does not indicate that the species is absent from the sampling point).
- _Medium:_ positive samples can interpreted as "Species DNA present in sample". There can be false negatives. 
- _High:_ positive samples can be interpreted as "Species DNA present in sample". There is low risk of false negatives.

**3. qPCR Experimental Methods:** You are then required to fill in simple information on the experimental methods used to run your qPCR plates. 

A standardised PDF report, providing a confidence assessment of your qPCR results, will be generated to send to other scientists and decision-makers. 

## Usage

**Step 1:** The tool requires a _very specific data format_ to function, so please download a modifiable data template from the home page (Figure 1). The template consists of a csv file with headers:

```
Plate	Well	Sample_Type	DNA_Sample	Replicate	Target_Cq	IPC_Cq	Std_Conc	Sampling_Point	Extraction_Batch	Volume_Water_Processed

```
```
Header	Details
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

**Step 2:** Fill out the Assay Validation table by carefully considering and answering each of the questions (1 = Yes, 0 = No; Figure 2).

**Step 3:** Fill out the qPCR Experimental Settings you used to generate your qPCR results (Figure 3).

**Step 4:** If you have correctly downloaded and used the data template to fill in your qPCR results (see Appendix 1 for an example), upload this in the 'Data Upload Section' (Figure 4). It is likely a 'warning' table will be generated. If so, you will need to do one of two things:
  - (A) An error appears: Something has gone wrong and you will not be able to generate a report (Figure 5). However, this could be something as simple as your columns headers no longer matching the data template. Check the error message and fix. 
  - (B) Warnings appear: One or more things have happened that you need to be aware of (Figure 6). You can still generate a report, however.

**Step 5:** If your dataset is compatable with the tool, a table and one or more figures (depending on the number of plates you used) will appear. These can be checked to ensure your data looks as expected. If your qPCR methods included using standard DNA samples of known concentrations, a standard curve plot will appear (Figure 7). If not, a plot showing the DNA samples and cycle threshold at which a signal was detected will be shown (Figure 8).  

**Step 6:** Enter a project title for the report. Additionally, the user may choose to select for the report to include an appendix with interpretation of the results at the qPCR replicate level, which will generate a table with the results inputted to the template file plus an additional column with an interpretation of each replicate result.

**Step 7:** Select generate report. The generated report will give the date and time the report was generated, followed by four sections: summary of results, assay validation and experimental methods, results and interpretations, and appendix. An example of a generated report can be seen in Appendix 2, but a brief description of each section is as follows:

1.	_Summary of results_ includes the assay validation level, if contamination was observed in any controls, if positive controls amplified as expected, and if standards were included in the uploaded data the associated R2 and PCR efficiency values.
2.	_Assay validation_ includes the table filled out by the user that was used to calculate the confidence in the assay (Low|Medium|High).
3.	_Experimental Methods_ includes a table with methods used to generate the qPCR results with defaults or changes made by the user.
4.	_Results and interpretations_ includes a figure, and tables summarising the results at 1) the Sampling_Point level 2) the DNA sample level and 3) [optional] the qPCR replicate level. Samples will either be classified as positive, negative, inconclusive, or tentative. A tentative is the result of weak signal (Cq value below the limit of detection).

The script is implemented as a Shiny app which has a publicly available dashboard [here](https://vidasolutions.shinyapps.io/TOOL_dashboard/)

## Common pitfalls to watch out for 
- **Header names:** It's very easy to change the names of the headers from what is in the data template, but these have to match. 
- **Number of plates used:** You must specify the number of plates that feature in your dataset.
- 

## **Disclaimers**
- **Work in progress:** COASTER remains a work in progress and we will fix new bugs if identified by users of the tool.
- **Interpretations:** ...
