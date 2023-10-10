# Dynamic-and-Concordance-assisted-Learning-for-Risk-Stratification

Executable codes of the manuscript:

"Dynamic and Concordance-assisted Learning for Risk Stratification with Application to Alzheimer's Disease" by WEN LI, RUOSHA LI, ZIDING FENG, JING NING

The following files are included:

1. main.R
	This R file contains the example code. It gives an example of applying the proposed method to a cohort of subjects with longitudinally measured markers. 

2. data.csv
	A similated longitudinal dataset that includes the following variables:
		id: subject id
		Y: observed survival time
		delta: censoring indicator
		X1, X2: two longitudinal markers
		t: measurement time

3. HelpFunc.R
	This R file contains functions that are needed in main.R

4. like_two_markers.cpp
	This cpp file contains functions that are needed in main.R
