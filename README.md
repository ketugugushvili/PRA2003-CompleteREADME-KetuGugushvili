# PRA2003 - Chemistry-Themed Project

**Week 7 deliverable:** Updated README documenting the complete PRA2003 project and the weekly assignments.  
**Student:** Ketu Gugushvili (i6357186)  
**Course:** Programming PRA2003 (MSP)  
**Project theme:** Chemistry

## Project overview

The project analyses simulations of molecular species produced during combustion at high temperatures using Python.
The process includes multiple steps: reading molecule records and calculating momentum, determining molecule frequencies and corresponding statistical uncertainties, and finally evaluating differences between molecule and isotope pairs, as well as how these differences vary with momentum.

## Input data and file structure

* Each input file represents repeated simulations of molecular species produced during combustion at high temperatures. 
* The header contains:
  - event ID: ID of the experiment or simulation run
  - number of molecules tracked: Total number of molecules observed in this run 
* Each row after the header gives:
  - 3D momentum components (px, py, pz in 10^(-23)kg m/s)
  - an integer ID representing the molecule or isotope 

```text
<event_id> <number_of_molecules>
<px> <py> <pz> <molecule_id>
<px> <py> <pz> <molecule_id>
...
<next_event_id> <number_of_molecules>
...
```

- `output-Set0.txt` is used by the Week 2 exercise.
- `output-Set1.txt` through `output-Set10.txt` are the ten subsamples used for the later analysis.

### Molecules analysed 
The table below displays molecule types that are going to be analysed and their respective IDs: 

| Molecule ID | Molecule type            |
|-------------|--------------------------|
| 211         | Carbon monoxide          |
| -211        | Carbon-13 monoxide       |
| 321         | Nitric oxide             |
| -321        | Ionised NO               |
| 2212        | Water                    | 
| -2212       | Heavy water (D2O)        |
| 3122        | Methane (CH4)            |
| -3122       | Methyl ion (CH3-)        |
| 3312        | Ethylene (C2H4)          |
| -3312       | Ionised ethylene (C2H3-) |
| 3334        | Ozone (O3)               |
| -3334       | Superoxide anion (O2-)   |


## Week 1 Deliverables: Research Questions  

Questions to be answered: 
1. **What are the average counts of each molecular species and their statistical uncertainties?**
2. **Is there any asymmetry between the normal and the variant molecule?**
   - Quantify!
4. **Is there any asymmetry as a function of their momentum?**
   - Quantify! 
   - Magnitude of total momentum calculated as: 
$$
p = \sqrt{p_x^2 + p_y^2 + p_z^2}.
$$

## Methods

The methods were developed across Weeks 2–5. Each stage builds on the previous one to answer the three research questions. The relevant Python script is identified below.

### Week 2 Deliverables: Data reading and momentum calculation

**File:** `Week_2(1).py`

The first programming exercise introduces reading a text file, extracting molecule IDs and the momentum components, converting strings into numerical values, and calculating the magnitude of momentum:

$$
p = \sqrt{p_x^2 + p_y^2 + p_z^2}.
$$

The script prints a table with the molecule ID and its calculated momentum. This establishes the data-reading and momentum calculations used again in the later analysis.

**Current limitation:** This introductory script skips the first header and expects all following lines to contain four molecule fields. It needs a more general event parser to work with files containing additional headers.

### Week 3 Deliverables: Average molecule frequency and Poisson uncertainty

**File:** `Week_3(2).py`

The program asks the user for an input filename, selects the twelve molecule IDs of interest, counts their occurrences, and calculates the average number of molecules per event and its Poisson statistical uncertainty. The results are printed in a table for all twelve species.

For molecule ID $m$:

$$
\bar N_m = \frac{N_m}{N_{\mathrm{events}}},
$$

where $N_m$ is the number of occurrences of that molecule and $N_{\mathrm{events}}$ is the number of events analysed.

Assuming Poisson counting statistics, the uncertainty on the average is:

$$
\sigma_{\bar N_m} = \frac{\sqrt{N_m}}{N_{\mathrm{events}}}.
$$

The script includes filename input, protection against empty input, and handling for missing files and file-access errors. **Current limitation:** The uploaded version counts events containing molecule records, rather than counting all event headers, so empty events may be missed. This should be corrected and checked against reference event totals.

### Week 4 Deliverables: Subsampling and overall molecule-pair asymmetry 

**File:** `Week_4(3).py` (also includes the Week 5 analysis)

The deliverable of Week 4 is to calculate the average number of each of the 12 molecule IDs per event. The entire sample of the 10 files should be analysed as sub-samples, separately, and be assigned statistical uncertainty from the spread of the results.  

The analysis also questions whether there is evidence of an asymmetry between each molecule pair by calculating their mean differences:
- 211 and -211
- 321 and -321
- 2212 and -2212
- 3122 and -3122
- 3312 and -3312
- 3334 and -3334

#### Combining the ten subsamples

The entire dataset was analysed as 10 separate sub-samples. 
For each sub-sample, the number of molecules associated with each of the 12 IDs was counted, and the average number of molecules per event was calculated. 

* The results from the sub-samples were then combined using a weighted average `np.average(values, weights=weights)`
  
$$
\bar N_m = \frac{\sum_{i=1}^{K} n_i\bar N_{m,i}}{\sum_{i=1}^{K}n_i}, \qquad K=10,
$$

* The events number($n_i$) in each sub-sample($i$) was used as the weight 
       - This was done to increase the contribution to the final average 
       of sub-samples with more events proportionally. 
       - In case the number of events is the same in all sub-samples, the weighted average is the same as the standard average. 


#### Statistical uncertainty

* To calculate statistical uncertainty, the sub-sampling method was used.

* For each sub-sample, the average number of molecules per event was calculated.
* Then, the **sample standard deviation** was calculated, representing the spread of these results with `np.std(values, ddof=1)`

$$
s_m = \sqrt{\frac{1}{K-1}\sum_{i=1}^{K}\left(\bar N_{m,i}-\langle\bar N_m\rangle\right)^2}.
$$

* This standard deviation was used as the statistical uncertainty of the analysis.

The full dataset is divided into ten subsamples. Each file is analysed separately to obtain an average frequency per event for each molecule. The results are then combined, and the spread between subsamples is used to describe statistical uncertainty.

#### Molecule-Pair Asymmetry Analysis

* To determine if a molecule was asymmetric to its ID pair or not, the difference between their averages was calculated for each sample:

$$
D_i = \bar N_{+,i}-\bar N_{-,i}.
$$

The overall difference is the event-weighted mean of these differences. Its uncertainty is the sample standard deviation of the ten paired differences:

$$
\sigma_D = SD(D_1,D_2,\ldots,D_K).
$$

* To further assess the importance of these differences, the absolute average difference was divided by its statistical uncertainty:
   - `|difference| / uncertainty` (in this file, this ratio will be referenced as significance)
$$
R=\frac{|\bar D|}{\sigma_D}.
$$

* This ratio represents how large the observed difference is compared to the statistical variation between the subsamples. 
* A larger ratio means that the observed difference is large compared with the subsample spread, therefore significant.
   - This ratio was used in this analysis as a descriptive measure to support the asymmetry comparison. 

### Week 5 Deliverables: Momentum-dependent asymmetry 

**File:** `Week_4(3).py`

The next question is whether the difference between each molecule and its variant changes with momentum. Instead of calculating only one difference for the complete momentum range, the program calculates a **separate average frequency difference for every momentum bin**.

The configured momentum range is 0–20, divided into 20 equal-width bins. The program calculates each molecule's momentum magnitude and counts molecules directly into the appropriate bin. This avoids storing every individual momentum value and reduces memory usage.

For subsample $i$ and momentum bin $j$, the average frequency per event is:

$$
f_{\pm,i,j}=\frac{N_{\pm,i,j}}{n_i}.
$$

The paired difference for that subsample and bin is:

$$
D_{i,j}=f_{+,i,j}-f_{-,i,j}.
$$

The graph's y-value is the mean paired difference across the ten subsamples:

$$
\bar D_j=\frac{1}{K}\sum_{i=1}^{K}D_{i,j},
$$

and the vertical error bar is the sample standard deviation:

$$
\sigma_j=SD(D_{1,j},D_{2,j},\ldots,D_{K,j}).
$$

**Reading the graphs:**

- **X-axis:** Momentum magnitude, divided into bins.
- **Y-axis:** Average frequency difference per event between the two molecules in each momentum bin.
- **Positive y-value:** The positive-ID molecule occurs more frequently in that bin, on average per event.
- **Negative y-value:** The negative-ID molecule occurs more frequently.
- **Zero:** Equal average frequencies for the pair in that bin.
- **Vertical error bars:** Subsample standard deviation of the paired differences.
- **Horizontal error bars:** Half the momentum-bin width, showing the bin interval rather than measurement uncertainty in momentum.

The momentum-bin mean is currently equally weighted across subsamples, while the overall Week 4 mean is event-weighted. The two are directly comparable when the subsamples have equal event counts; otherwise, this difference in weighting should be considered.

The script is designed to generate one plot per pair, named `momentum_asymmetry_<positiveID>_<negativeID>.png`. Once the corrected code has been run and checked, the plots can be added here. A single unusual bin among twenty does not automatically establish a momentum-dependent effect.



## Average Number of Molecules per Event

The results obtained from the complete sample are:

-------------*----------------------*---------------------------
Molecule ID  | Average per event    | Statistical uncertainty
-------------*----------------------*---------------------------
211          | 19.94950794          | 0.03274619
-211         | 19.91721122          | 0.03187334
321          | 2.50914781           | 0.00477369
-321         | 2.50345680           | 0.00550369
2212         | 1.20803410           | 0.00189633
-2212        | 1.18416095           | 0.00241216
3122         | 0.27659851           | 0.00107440
-3122        | 0.27169595           | 0.00098538
3312         | 0.03944116           | 0.00028368
-3312        | 0.03900049           | 0.00040248
3334         | 0.00118710           | 0.00004168
-3334        | 0.00115158           | 0.00005079


--------------*----------------------*--------------------*---------------
 Asymmetry Between Molecule Pairs 
--------------*----------------------*--------------------*---------------
  Pair        |  Average difference  |  Uncertainty       |  Significance  
--------------*----------------------*--------------------*---------------
 211 vs -211  |           0.03229671 |         0.00451870 |         7.15
 321 vs -321  |           0.00569100 |         0.00329367 |         1.73
2212 vs -2212 |           0.02387314 |         0.00237277 |         10.06
3122 vs -3122 |           0.00490256 |         0.00058316 |         8.41
3312 vs -3312 |           0.00044067 |         0.00048793 |         0.90
3334 vs -3334 |           0.00003551 |         0.00006854 |         0.52


## Interpretation of the Molecule Pairs

### 211 vs -211

The average difference is: `0.03229671 ± 0.00451870`

The difference is approximately 7.15 times the statistical uncertainty. 
--> Thus, the observed difference is significant due to how large it 
is relative to the statistical variation between the sub-samples.
----> The results provide evidence of an asymmetry between 211 and -211.


### 321 vs -321

The average difference is: `0.00569100 ± 0.00329367`

The difference is approximately 1.73 times the statistical uncertainty.
--> Thus, the difference is not significantly large compared with the 
statistical variation between the sub-samples. 
----> The evidence doesn't support an asymmetry between 321 and -321 completely.


### 2212 vs -2212

The average difference is: `0.02387314 ± 0.00237277`

The difference is approximately 10.06 times the statistical uncertainty.
--> The observed difference is much larger than the statistical variation. 
----> The results provide evidence of an asymmetry between 2212 and -2212.


### 3122 vs -3122

The average difference is: `0.00490256 ± 0.00058316`

The difference is approximately 8.41 times the statistical uncertainty.
--> The difference is significantly large compared with the statistical variation.
----> The results provide evidence of an asymmetry between 3122 and -3122.


### 3312 vs -3312

The average difference is: `0.00044067 ± 0.00048793`

The difference is approximately 0.90 times the statistical uncertainty.
--> The difference is way too small relative to the statistical variation.   
----> There is no evidence of an asymmetry between 3312 and -3312.


### 3334 vs -3334

The average difference is: `0.00003551 ± 0.00006854`

The difference is approximately 0.52 times the statistical uncertainty.
--> The difference is smaller than its uncertainty.
----> There is no evidence of an asymmetry between 3334 and -3334.

 
