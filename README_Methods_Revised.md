# PRA2003 — Chemistry-Themed Programming Project

**Week 7 deliverable:** Updated README documenting the complete PRA2003 project and its development across the weekly assignments.  
**Student:** Ketu Gugushvili (i6357186)  
**Course:** Programming PRA2003 (MSP)  
**Project theme:** Chemistry

## Project overview

This project uses Python to analyse simulated molecular species produced during combustion at high temperatures. The work develops step by step: first reading molecule records and calculating momentum, then determining molecule frequencies and statistical uncertainties, and finally investigating differences between related molecule pairs, including how these differences vary with momentum.

The dataset consists of simulation **events**, each containing a header followed by molecule records. Each molecule is identified by an integer ID and three momentum components. The analyses focus on twelve IDs, organised into six pairs.

### Molecules analysed

| Molecule ID | Molecule type |
|---:|---|
| 211 | Carbon monoxide |
| -211 | Carbon-13 monoxide |
| 321 | Nitric oxide |
| -321 | Ionised NO |
| 2212 | Water |
| -2212 | Heavy water (D2O) |
| 3122 | Methane (CH4) |
| -3122 | Methyl ion (CH3-) |
| 3312 | Ethylene (C2H4) |
| -3312 | Ionised ethylene (C2H3-) |
| 3334 | Ozone (O3) |
| -3334 | Superoxide anion (O2-) |

*These names and ID assignments follow the project's supplied mapping; the IDs should not be interpreted as standard chemical identifiers outside this dataset.*

## Research questions (Week 1)

The project aims to answer three questions:

1. **What are the average counts of each molecular species and their statistical uncertainties?**
2. **Is there any asymmetry between a normal molecule and its variant?** Quantify the difference and its uncertainty.
3. **Is there any asymmetry as a function of momentum?** Quantify the average frequency difference for each momentum interval and its uncertainty.

The third question uses the magnitude of the three-dimensional momentum vector:

$$
p = \sqrt{p_x^2 + p_y^2 + p_z^2}.
$$

## Input data and file structure

The input data are text files containing repeated simulation events. Each event has a two-integer header: the event ID and the number of molecule records associated with that event. The following rows contain three momentum components and one molecule ID.

```text
<event_id> <number_of_molecules>
<px> <py> <pz> <molecule_id>
<px> <py> <pz> <molecule_id>
...
<next_event_id> <number_of_molecules>
...
```

The supplied project description gives momentum components in units of $10^{-23}\,\mathrm{kg\,m/s}$. The axis labels in the current Week 4/5 script use `1e-20 kg m/s`; **confirm the intended scale with the dataset documentation and make the labels consistent before publication**. The code itself calculates the vector magnitude without a unit conversion.

- `output-Set0.txt` is used by the Week 2 exercise.
- `output-Set1.txt` through `output-Set10.txt` are the ten subsamples used for the later analysis.
- A file may contain many events; a file is **not** the same as a single event.

## Methods

The methods were developed across Weeks 2–5. Each stage builds on the previous one to answer the three research questions. The relevant Python script is identified below.

### Data reading and momentum calculation (Week 2)

**File:** `Week_2(1).py`

The first programming exercise introduces reading a text file, extracting molecule IDs and the momentum components, converting strings into numerical values, and calculating the magnitude of momentum:

$$
p = \sqrt{p_x^2 + p_y^2 + p_z^2}.
$$

The script prints a table with the molecule ID and its calculated momentum. This establishes the data-reading and momentum calculations used again in the later analysis.

**Current limitation:** This introductory script skips the first header and expects all following lines to contain four molecule fields. It needs a more general event parser to work with files containing additional headers.

### Average molecule frequency and Poisson uncertainty (Week 3)

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

### Subsampling and overall molecule-pair asymmetry (Week 4)

**File:** `Week_4(3).py` (also includes the Week 5 analysis)

The full dataset is divided into ten subsamples. Each file is analysed separately to obtain an average frequency per event for each molecule. The results are then combined, and the spread between subsamples is used to describe statistical uncertainty.

#### Combining the ten subsamples

The combined average is weighted by the number of events in each subsample:

$$
\bar N_m = \frac{\sum_{i=1}^{K} n_i\bar N_{m,i}}{\sum_{i=1}^{K}n_i}, \qquad K=10,
$$

where $n_i$ is the event count in subsample $i$. This gives larger subsamples proportionally more weight. If all subsamples have equal event counts, the weighted and unweighted means agree.

The reported subsampling uncertainty is the **sample standard deviation** of the ten subsample averages, calculated with `np.std(values, ddof=1)`:

$$
s_m = \sqrt{\frac{1}{K-1}\sum_{i=1}^{K}\left(\bar N_{m,i}-\langle\bar N_m\rangle\right)^2}.
$$

This measures the spread between subsample results; it is **not** the standard error of the combined mean.

#### Overall molecule-pair asymmetry

The six comparisons are `211/-211`, `321/-321`, `2212/-2212`, `3122/-3122`, `3312/-3312`, and `3334/-3334`.

For each pair, the difference is calculated **within each subsample**:

$$
D_i = \bar N_{+,i}-\bar N_{-,i}.
$$

The overall difference is the event-weighted mean of these differences. Its uncertainty is the sample standard deviation of the ten paired differences:

$$
\sigma_D = SD(D_1,D_2,\ldots,D_K).
$$

Calculating paired differences first allows variation shared by both molecule types within a subsample to be reflected in the uncertainty. Simply adding the separate molecule uncertainties in quadrature would assume zero covariance and is not the method used here.

To compare the size of each observed difference with the subsample spread, the project also calculates:

$$
R=\frac{|\bar D|}{\sigma_D}.
$$

In the program this is labelled *significance*, but here it is used as a **descriptive difference-to-spread ratio**, not as a formal Gaussian significance or hypothesis-test result.

### Momentum-dependent asymmetry (Week 5)

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

The momentum-bin mean is currently equally weighted across subsamples, while the overall Week 4 mean is event-weighted. The two are directly comparable when the subsamples have equal event counts; otherwise this difference in weighting should be considered.

The script is designed to generate one plot per pair, named `momentum_asymmetry_<positiveID>_<negativeID>.png`. Once the corrected code has been run and checked, the plots can be added here. A single unusual bin among twenty does not automatically establish a momentum-dependent effect.

## Results and interpretation

### Previously obtained overall results

These results were recorded during the earlier analysis. They should be checked again after the latest event-parsing changes.

| Molecule ID | Average per event | Subsampling uncertainty |
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

| Molecule pair | Average difference per event | Uncertainty | Difference / uncertainty |
|---|---:|---:|---:|
| 211 vs -211 | 0.03229671 | 0.00451870 | 7.15 |
| 321 vs -321 | 0.00569100 | 0.00329367 | 1.73 |
| 2212 vs -2212 | 0.02387314 | 0.00237277 | 10.06 |
| 3122 vs -3122 | 0.00490256 | 0.00058316 | 8.41 |
| 3312 vs -3312 | 0.00044067 | 0.00048793 | 0.90 |
| 3334 vs -3334 | 0.00003551 | 0.00006854 | 0.52 |

### Interpretation of the overall asymmetry

The differences for **211/-211**, **2212/-2212**, and **3122/-3122** are relatively large compared with the observed subsample spread (ratios 7.15, 10.06, and 8.41). These pairs show the strongest descriptive asymmetries in the reported results.

For **321/-321**, the difference is 1.73 times the reported uncertainty, providing weaker evidence of a consistent difference. For **3312/-3312** and **3334/-3334**, the ratios are 0.90 and 0.52, respectively: their observed differences are smaller than the reported spread and are consistent with approximate symmetry at this descriptive level.

These comparisons do not, by themselves, establish formal statistical significance. A formal test would require an appropriate sampling-error model and consideration of multiple comparisons.

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

Clone or download the repository and open a terminal in its folder. Optionally create and activate a virtual environment:

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

The filenames, momentum range and number of bins are configured near the top of `Week_4(3).py` using `FILENAMES`, `MOMENTUM_MIN`, `MOMENTUM_MAX` and `NUMBER_OF_BINS`.

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

## Contributors

Ketu Gugushvili — PRA2003 project README and programming work.  
Python 2 project group — add other contributors' names or GitHub handles as appropriate.

## Acknowledgements and references

- PRA2003 programming course materials and supplied chemistry-themed simulation dataset.
- [Python documentation](https://docs.python.org/3/)
- [NumPy documentation](https://numpy.org/doc/)
- [Matplotlib documentation](https://matplotlib.org/stable/)

## License and data availability

No open-source license has been specified for this repository. The course-provided datasets may have separate redistribution restrictions and should not be published without permission.
