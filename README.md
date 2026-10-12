# PRA2003 - Chemistry-Themed Project

**Week 7 deliverable:** Updated README documenting the complete PRA2003 project and the weekly assignments.  
**Student:** Ketu Gugushvili (i6357186)  
**Course:** Programming PRA2003 (MSP)  
**Project theme:** Chemistry

## Project overview

* The project analyses simulations of molecular species produced during combustion at high temperatures using Python.
* The process includes multiple steps: reading molecule records and calculating momentum, determining molecule frequencies and corresponding statistical uncertainties, and finally evaluating differences between molecule and isotope pairs, as well as how these differences vary with momentum.

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

### Week 2: Data reading and momentum calculation

**File:** `Week_2.py`

The Week 2 code reads a text file, selects molecule IDs and the momentum components, converts strings into numerical values, and calculates the magnitude of momentum:

$$
p = \sqrt{p_x^2 + p_y^2 + p_z^2}.
$$

The results are displayed in a table of the molecule ID and its calculated momentum.

### Week 3: Average molecule frequency and Poisson uncertainty

**File:** `Week_3.py`

The Week 3 code asks the user for an input filename, selects the 12 molecule IDs of interest, counts their occurrences, and calculates the average number of molecules per event and its Poisson statistical uncertainty. 

For molecule ID $m$:

$$
\bar N_m = \frac{N_m}{N_{\mathrm{events}}},
$$

where $N_m$ is the number of occurrences of that molecule and $N_{\mathrm{events}}$ is the number of events analysed.

Assuming Poisson counting statistics, the uncertainty on the average is:

$$
\sigma_{\bar N_m} = \frac{\sqrt{N_m}}{N_{\mathrm{events}}}.
$$

The results for each of the 12 molecule IDs are displayed in the table. 

### Week 4: Subsampling and overall molecule-pair asymmetry 

**File:** `Week_4(3).py` (also includes the Week 5 analysis)

The Week 4 code calculates the average number of each of the 12 molecule IDs per event. The entire sample of the 10 files is analysed as sub-samples, separately, and assigned statistical uncertainty from the spread of the results.  

The analysis also questions whether there is evidence of an asymmetry between each molecule pair by calculating their mean differences:
- 211 and -211
- 321 and -321
- 2212 and -2212
- 3122 and -3122
- 3312 and -3312
- 3334 and -3334

#### Combining the ten subsamples
 
* For each sub-sample, the code counts the number of molecules associated with each of the 12 IDs and calculates the average number of molecules per event. 
* Then, it combines the results from the subsamples using a weighted average `np.average(values, weights=weights)`:
  
$$
\bar N_m = \frac{\sum_{i=1}^{K} n_i\bar N_{m,i}}{\sum_{i=1}^{K}n_i}, \qquad K=10,
$$

* The events number($n_i$) in each sub-sample($i$) is used as the weight
  - This is done to increase the contribution to the final average  of sub-samples with more events proportionally.
  - In case the number of events is the same in all sub-samples, the weighted average is the same as the standard average. 

#### Statistical uncertainty

* The program uses the sub-sampling method to calculate statistical uncertainty.
* For each sub-sample, the average number of molecules per event is calculated.
* Then, statistical uncertainty is calculated as the **sample standard deviation**, representing the spread of these results with `np.std(values, ddof=1)`:
  
$$
s_m = \sqrt{\frac{1}{K-1}\sum_{i=1}^{K}\left(\bar N_{m,i}-\langle\bar N_m\rangle\right)^2}.
$$

#### Molecule-Pair Asymmetry Analysis

* To determine if a molecule was asymmetric to its ID pair or not, the difference between their averages is calculated for each sample:
  - Since data for each ID comes from the same dataset, it isn't possible to assume no correlation between them
  - The difference that doesn't assume no correlation is expressed as:
    
$$
D_i = \bar N_{+,i}-\bar N_{-,i}.
$$

* The overall difference is the event-weighted mean of these differences.
* Its statistical uncertainty is the sample standard deviation of the ten paired differences:

$$
SD = SD(D_1,D_2,\ldots,D_K).
$$

* To further assess the importance of these differences, the absolute average difference is divided by its statistical uncertainty:
   - `|difference| / uncertainty` (in this file, this ratio will be referenced as significance) 

$$
\sigma =\frac{|\bar D|}{SD}.
$$

* This ratio represents how large the observed difference is compared to the statistical variation between the subsamples. 
* A larger ratio means that the observed difference is large compared with the subsample spread, therefore significant.
   - This ratio was used in this analysis as a descriptive measure to support the asymmetry comparison:
     - 1$\sigma$: 1 in 6 chance the result is caused by fluctuations
     - 2$\sigma$: 1 in 44 chance
     - 3$\sigma$: correspond to a 1-in-740 chance of a statistical quirk
     - 4$\sigma$: corresponds to a 1-in-32,000 chance
     - 5$\sigma$: is a 1-in-3.5 million chance

### Week 5 Deliverables: Momentum-dependent asymmetry 

**File:** `Week_4.py`

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
SD=SD(D_{1,j},D_{2,j},\ldots,D_{K,j}).
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

## Results and interpretation

### Week 4: Average Number of Molecules per Event

The results obtained from the complete sample are:

| Molecule ID | Average per event | Statistical uncertainty |
|---:|---:|---:|
| 211 | 19.94950794 | 0.03274619 |
| -211 | 19.91721122 | 0.03187334 |
| 321 | 2.50914781 | 0.00477369 |
| -321 | 2.50345680 | 0.00550369 |
| 2212 | 1.20803410 | 0.00189633 |
| -2212 | 1.18416095 | 0.00241216 |
| 3122 | 0.27659851 | 0.00107440 |
| -3122 | 0.27169595 | 0.00098538 |
| 3312 | 0.03944116 | 0.00028368 |
| -3312 | 0.03900049 | 0.00040248 |
| 3334 | 0.00118710 | 0.00004168 |
| -3334 | 0.00115158 | 0.00005079 |

| Molecule pair | Average difference per event | Statistical Uncertainty | Significance | Asymmetry |
|---|---:|---:|---:|---:|
| 211 vs -211 | 0.03229671 | 0.00451870 | 7.15 | Yes |
| 321 vs -321 | 0.00569100 | 0.00329367 | 1.73 | No |
| 2212 vs -2212 | 0.02387314 | 0.00237277 | 10.06 | Yes |
| 3122 vs -3122 | 0.00490256 | 0.00058316 | 8.41 | Yes |
| 3312 vs -3312 | 0.00044067 | 0.00048793 | 0.90 | No |
| 3334 vs -3334 | 0.00003551 | 0.00006854 | 0.52 | No |


#### Interpretation of the Molecule Pairs

##### 211 vs -211

The average difference is: `0.03229671 ± 0.00451870`

* The difference is approximately 7.15 times the statistical uncertainty.
  - Thus, the observed difference is significant due to how large it is relative to the statistical variation between the sub-samples.
  - The results provide evidence of an asymmetry between 211 and -211.

##### 321 vs -321

The average difference is: `0.00569100 ± 0.00329367`

* The difference is approximately 1.73 times the statistical uncertainty.
  - Thus, the difference is not significantly large compared with the statistical variation between the sub-samples.
  - The evidence doesn't support an asymmetry between 321 and -321 completely.

##### 2212 vs -2212

The average difference is: `0.02387314 ± 0.00237277`

* The difference is approximately 10.06 times the statistical uncertainty.
  - The observed difference is much larger than the statistical variation.
  - The results provide evidence of an asymmetry between 2212 and -2212.

##### 3122 vs -3122

The average difference is: `0.00490256 ± 0.00058316`

* The difference is approximately 8.41 times the statistical uncertainty.
  - The difference is significantly large compared with the statistical variation.
  - The results provide evidence of an asymmetry between 3122 and -3122.

##### 3312 vs -3312

The average difference is: `0.00044067 ± 0.00048793`

* The difference is approximately 0.90 times the statistical uncertainty.
  - The difference is way too small relative to the statistical variation.
  - There is no evidence of an asymmetry between 3312 and -3312.

##### 3334 vs -3334

The average difference is: `0.00003551 ± 0.00006854`

* The difference is approximately 0.52 times the statistical uncertainty.
  - The difference is smaller than its uncertainty.
  - There is no evidence of an asymmetry between 3334 and -3334.


### Week 5: Momentum-Dependent Asymmetry

#### 211 vs -211: Carbon monoxide and Carbon-13 monoxide

<img src="Plots/momentum_asymmetry_211_-211.png" alt="Description of image" width="600">

#### 321 vs -321: Nitric Oxide and ionised NO

<img src="Plots/momentum_asymmetry_321_-321.png" alt="Description of image" width="600">

#### 2212 vs -2212: Water and heavy water 

<img src="Plots/momentum_asymmetry_2212_-2212.png" alt="Description of image" width="600">

#### 3122 vs -3122: Methane and methyl ion

<img src="Plots/momentum_asymmetry_3122_-3122.png" alt="Description of image" width="600">

#### 3312 vs -3312: Ethylene and ionised ethylene

<img src="Plots/momentum_asymmetry_3312_-3312.png" alt="Description of image" width="600">

#### 3334 vs -3334: Ozone  and superoxide anion 

<img src="Plots/momentum_asymmetry_3334_-3334.png" alt="Description of image" width="600">

##### Overall comparison and interpretation of momentum-dependent asymmetry 

| Molecule pair | Pattern visible in momentum plot | Interpretation |
|---|---|---|
| 211 vs -211 | Positive difference, largest at low momentum, decreasing with momentum | Clear descriptive low-momentum concentration |
| 321 vs -321 | Small differences, variable sign at low momentum; uncertainties often span zero | No clear consistent momentum trend |
| 2212 vs -2212 | Positive difference that generally decreases with momentum | Clear descriptive low-momentum concentration |
| 3122 vs -3122 | Mostly positive, larger at low momentum, but substantial spread | Suggestive decreasing pattern with considerable uncertainty |
| 3312 vs -3312 | Very small differences, mostly consistent with zero within plotted spread | No clear momentum-dependent asymmetry |
| 3334 vs -3334 | Tiny differences of mixed sign with large relative spread | No clear momentum-dependent asymmetry |


## Repository structure

```text
PRA2003-project/
├── README.md
├── Week_2(1).py                 # Week 2: momentum calculation
├── Week_3(2).py                 # Week 3: one-file molecule frequencies
├── Week_4(3).py                 # Weeks 4–5: subsamples, pair and momentum analysis
├── output-Set0.txt              # Week 2 input, if available
├── output-Set1.txt              # Subsample 1
├── ...
├── output-Set10.txt             # Subsample 10
└── momentum_asymmetry_*.png     # Generated momentum plots
```

Large course-provided data files may be excluded from the public repository. In that case, obtain them through the course and place them in the working directory before running the scripts.

## Requirements and installation

- Python 3.9 or newer
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- Python standard-library modules `os` and `math` (no separate installation required)

Clone or download the repository and open a terminal in its folder. Optionally, create and activate a virtual environment:

```bash
python -m venv .venv
```

**macOS / Linux:**

```bash
source .venv/bin/activate
```

**Windows PowerShell:**

```powershell
.venv\Scripts\Activate.ps1
```

Install the external packages:

```bash
python -m pip install numpy matplotlib
```

Place the required `.txt` files alongside the scripts, or adjust the filenames in the code to point to their location.

## How to run the project

Run the scripts from the project directory in the order corresponding to the weekly work.

**Week 2 — Calculate molecule momentum:**

```bash
python "Week_2(1).py"
```

The script reads `output-Set0.txt` and prints molecule IDs and momentum magnitudes.

**Week 3 — Analyse one file:**

```bash
python "Week_3(2).py"
```

Enter the input filename when prompted, for example `output-Set1.txt`. The script prints the number of detected events and the average frequency and Poisson uncertainty for each of the twelve molecule IDs.

**Weeks 4–5 — Analyse all ten subsamples and generate plots:**

```bash
python "Week_4(3).py"
```

This script expects `output-Set1.txt` through `output-Set10.txt`. It prints the per-subsample and combined frequency tables, overall molecule-pair differences, and momentum-dependent bin results, then saves six PNG graphs.

The filenames, momentum range, and number of bins are configured near the top of `Week_4(3).py` using `FILENAMES`, `MOMENTUM_MIN`, `MOMENTUM_MAX` and `NUMBER_OF_BINS`.

**Status of the uploaded code:** Before relying on the output, correct the event-header parsing in `Week_4(3).py` and the indentation of the momentum-bin results printing loop. The current Week 3 event counting and Week 2 multi-event parsing also require review. The scripts have not been verified against the complete input dataset in this README.

## Program output

The project produces:

1. A Week 2 console table of molecule IDs and momentum magnitudes.
2. A Week 3 table of average counts per event and Poisson uncertainties for one input file.
3. Weeks 4–5 tables of individual-subsample results, combined averages, pair differences and momentum-bin differences.
4. Six momentum-dependent PNG graphs (one for each molecule pair).
5. Diagnostics such as the number of analysed events and the count of selected molecules outside the chosen momentum range.

The current scripts do not export numerical tables to CSV.

## Validation, assumptions and limitations

- Event headers must be counted correctly, including events with zero molecule records. Compare the total number of events with the course-provided reference.
- Malformed input records and non-integer molecule IDs should be rejected rather than silently interpreted as valid data.
- Week 2's simple reader is intended as an introductory exercise and is not yet a general multi-event parser.
- The Poisson uncertainty in Week 3 assumes Poisson counting statistics.
- The subsample standard deviation in Weeks 4–5 measures spread between subsamples; it is not automatically the standard error of the combined estimate.
- Pairwise differences are calculated within the same subsample to retain the effect of shared variation.
- Molecules outside the configured momentum range are excluded from the momentum-bin plots; the diagnostic count should be checked.
- Momentum units and graph-axis labels must be made consistent with the supplied dataset.
- The numerical results shown above were obtained earlier and must be rechecked after the latest parsing corrections.

## Reproducibility and future improvements

To reproduce the analysis, use the same input files, filenames, Python packages and momentum-bin settings. Check the event counts first, then compare the individual molecule averages and pair differences with the previously recorded results before interpreting the graphs.

Potential improvements include a shared event parser across all three scripts, automated tests for empty events and malformed records, CSV export for numerical tables, configurable input paths and verified momentum plots embedded in this README.

## References

- PRA2003 programming course materials and supplied chemistry-themed simulation dataset.
- [Python documentation](https://docs.python.org/3/)
- [NumPy documentation](https://numpy.org/doc/)
- [Matplotlib documentation](https://matplotlib.org/stable/)


