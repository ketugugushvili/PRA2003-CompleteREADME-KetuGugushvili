############################################################
# Week 4 & 5: Complete molecule project analysis 
############################################################

import os
import numpy as np
import matplotlib.pyplot as plt
import csv

############################################################
# The parameters for the analysis
############################################################

MOLECULE_PAIRS = [
    (211, -211),
    (321, -321),
    (2212, -2212),
    (3122, -3122),
    (3312, -3312),
    (3334, -3334)
]

MOLECULES_OF_INTEREST = [
    molecule_id
    for pair in MOLECULE_PAIRS
    for molecule_id in pair
]


# Automatically generate filenames for Sets 1-10.
FILENAMES = [
    f"output-Set{i}.txt"
    for i in range(1, 11)
]


# Momentum analysis settings.
MOMENTUM_MIN = 0
MOMENTUM_MAX = 20
NUMBER_OF_BINS = 20

if MOMENTUM_MAX <= MOMENTUM_MIN:
    raise ValueError("Invalid momentum range.")

if NUMBER_OF_BINS <= 0:
    raise ValueError("Number of bins must be positive.")


############################################################
# Check and remove any empty data lines without molecules
############################################################

def is_molecule_line(tokens):

    if len(tokens) != 4:
        return False

    try:
        px, py, pz = map(float, tokens[:3])
        molecule_id = float(tokens[3])

        if not np.all(np.isfinite([px, py, pz, molecule_id])):
            return False

        if not molecule_id.is_integer():
            return False

    except ValueError:
        return False

    return True


############################################################
# Analyse one sub-sample
############################################################

def analyse_subsample(filename):

    # Initialise molecule counts.
    molecule_counts = {
        molecule_id: 0
        for molecule_id in MOLECULES_OF_INTEREST
    }

    # Initialise momentum-bin counts for each molecule.
    momentum_counts = {
        molecule_id: np.zeros(NUMBER_OF_BINS, dtype=int)
        for molecule_id in MOLECULES_OF_INTEREST
    }

    bin_edges = create_momentum_bins()

    number_of_events = 0
    outside_range = 0
    event_started = False

    with open(filename, "r") as data_file:

        for line_number, line in enumerate(data_file, start=1):

            tokens = line.split()

            if not tokens:
                continue


            if len(tokens) == 2:

                # Count the previous event only if
                # it contained molecule records.
                if event_started:
                    number_of_events += 1

                event_started = False

                # Validate the new event header.
                try:
                    int(tokens[0])
                    int(tokens[1])
                except ValueError:
                    raise ValueError(
                        f"Invalid event header on line {line_number}"
                    )

            elif is_molecule_line(tokens):

                px, py, pz = map(float, tokens[:3])
                molecule_id = int(float(tokens[3]))

                # Mark the event as nonempty.
                event_started = True

                if molecule_id in molecule_counts:

                    molecule_counts[molecule_id] += 1

                    momentum = np.sqrt(
                        px**2 + py**2 + pz**2
                    )

                    bin_index = np.searchsorted(
                        bin_edges,
                        momentum,
                        side="right"
                    ) - 1

                    if momentum == MOMENTUM_MAX:
                        bin_index = NUMBER_OF_BINS - 1

                    if 0 <= bin_index < NUMBER_OF_BINS:
                        momentum_counts[molecule_id][bin_index] += 1
                    else:
                        outside_range += 1

            else:
                raise ValueError(
                    f"Invalid data on line {line_number}: {line.strip()}"
                )


    # Count the final event.

    if event_started:
        number_of_events += 1

    if number_of_events == 0:
        raise ValueError(
            f"No events were found in {filename}."
        )

    averages = {}
    poisson_uncertainties = {}

    for molecule_id in MOLECULES_OF_INTEREST:

        total_count = molecule_counts[molecule_id]

        averages[molecule_id] = (
            total_count / number_of_events
        )

        poisson_uncertainties[molecule_id] = (
            np.sqrt(total_count) / number_of_events
        )

    print(
        "Molecules outside momentum range:",
        outside_range
    )

    return (
        averages,
        poisson_uncertainties,
        momentum_counts,
        number_of_events
    )

############################################################
# Calculate the weighted average
############################################################
def weighted_average(values, weights):
    # The number of events in each sub-sample is used
    # as its weight.
    average = np.average(
    values,
    weights=weights
    )
    return average

############################################################
# Calculate the sub-sampling uncertainty
############################################################

def calculate_subsample_uncertainty(values):

    # At least two sub-samples are needed
    # to calculate a sample standard deviation.
    if len(values) < 2:
        return 0.0


    return np.std(
        values,
        ddof=1
    )


############################################################
# Calculate significance of the measurements
############################################################

def calculate_significance(
    difference,
    uncertainty
):

    # Protection against division by zero.
    if uncertainty == 0:

        if difference == 0:
            return 0.0

        return np.inf


    return (
        abs(difference)
        / uncertainty
    )


############################################################
# Analyse all ten sub-samples
############################################################

def analyse_all_subsamples(filenames):

    # Store average for every ID from every file.
    subsample_averages = {
        molecule_id: []
        for molecule_id in MOLECULES_OF_INTEREST
    }


    # Store Poisson uncertainty for every ID
    # from every individual file.
    subsample_uncertainties = {
        molecule_id: []
        for molecule_id in MOLECULES_OF_INTEREST
    }


    # Momentum information from every file.
    all_momentum_counts = [] 


    # Number of events in every sub-sample.
    events_per_subsample = []


    for filename in filenames:

        print()
        print(
            "Analysing",
            filename,
            "..."
        )


        (
            averages,
            poisson_uncertainties,
            momentum_counts,
            number_of_events
        ) = analyse_subsample(filename)


        events_per_subsample.append(
            number_of_events
        )


        all_momentum_counts.append(momentum_counts)


        for molecule_id in MOLECULES_OF_INTEREST:

            subsample_averages[
                molecule_id
            ].append(
                averages[molecule_id]
            )


            subsample_uncertainties[
                molecule_id
            ].append(
                poisson_uncertainties[
                    molecule_id
                ]
            )


        print(
            "Number of events:",
            number_of_events
        )


    return (
        subsample_averages,
        subsample_uncertainties,
        all_momentum_counts,
        events_per_subsample
    )


############################################################
# Calculate the final averages 
############################################################

def calculate_final_results(
    subsample_averages,
    events_per_subsample
):

    final_results = {}


    for molecule_id in MOLECULES_OF_INTEREST:

        values = np.array(
            subsample_averages[
                molecule_id
            ]
        )


        # Weighted average of the sub-samples.
        final_average = np.average(
            values,
            weights=events_per_subsample
        )


        # Statistical uncertainty from the
        # spread of the sub-sample results.
        uncertainty = (
            calculate_subsample_uncertainty(
                values
            )
        )


        final_results[
            molecule_id
        ] = (
            final_average,
            uncertainty
        )


    return final_results


############################################################
# Calculate overall molecule pair asymmetries 
############################################################


def analyse_pairs(subsample_averages, events_per_subsample):

    pair_results = {}

    for positive_id, negative_id in MOLECULE_PAIRS:

        differences = []

        for i in range(len(events_per_subsample)):

            positive_average = subsample_averages[positive_id][i]
            negative_average = subsample_averages[negative_id][i]

            differences.append(
                positive_average - negative_average
            )

        # Same weighted average as Week_4(4).py
        average_difference = weighted_average(
            differences,
            events_per_subsample
        )

        # Same subsampling uncertainty
        difference_uncertainty = calculate_subsample_uncertainty(
            differences
        )

        significance = calculate_significance(
            average_difference,
            difference_uncertainty
        )

        pair_results[(positive_id, negative_id)] = (
            average_difference,
            difference_uncertainty,
            significance
        )

    return pair_results


############################################################
# Create momentum bins
############################################################

def create_momentum_bins():

    return np.linspace(
        MOMENTUM_MIN,
        MOMENTUM_MAX,
        NUMBER_OF_BINS + 1
    )


############################################################
# Analyse momentum-dependent average frequency difference
############################################################

def analyse_momentum_asymmetry(
    all_momentum_counts,
    events_per_subsample
):

    bin_edges = create_momentum_bins()
    momentum_results = {}

    if len(all_momentum_counts) != len(events_per_subsample):
        raise ValueError("Subsample data and event counts do not match.")

    if len(events_per_subsample) < 2:
        raise ValueError("At least two subsamples are required.")

    for positive_id, negative_id in MOLECULE_PAIRS:

        subsample_differences = []

        for counts, number_of_events in zip(
            all_momentum_counts,
            events_per_subsample
        ):

            if number_of_events <= 0:
                raise ValueError("Invalid number of events.")

            # Frequency per event for each momentum bin.
            positive_frequency = (
                counts[positive_id] / number_of_events
            )

            negative_frequency = (
                counts[negative_id] / number_of_events
            )

            # Paired difference within this subsample.
            difference = (
                positive_frequency - negative_frequency
            )

            subsample_differences.append(difference)

        subsample_differences = np.array(
            subsample_differences
        )

        # Mean difference for each momentum bin.
        mean_differences = np.mean(
            subsample_differences,
            axis=0
        )

        # Standard deviation of paired differences.
        uncertainties = np.std(
            subsample_differences,
            axis=0,
            ddof=1
        )

        significances = np.array([
            calculate_significance(difference, uncertainty)
            for difference, uncertainty in zip(
                mean_differences,
                uncertainties
            )
        ])

        momentum_results[
            (positive_id, negative_id)
        ] = (
            mean_differences,
            uncertainties,
            significances
        )

    return momentum_results, bin_edges


############################################################
# Print results for each individual sub-sample
############################################################

def print_subsample_results(
    subsample_averages,
    subsample_uncertainties
):

    print()
    print("=" * 82)
    print("RESULTS FOR EACH INDIVIDUAL SUB-SAMPLE")
    print("=" * 82)

    # Table headings
    print(
        f"{'Sub-sample':>12} | "
        f"{'Molecule ID':>12} | "
        f"{'Average per event':>20} | "
        f"{'Statistical uncertainty':>24}"
    )

    print("-" * 82)


    # Find how many sub-samples were analysed.
    number_of_subsamples = len(
        subsample_averages[
            MOLECULES_OF_INTEREST[0]
        ]
    )


    # Print one section for each sub-sample.
    for i in range(number_of_subsamples):

        # Create the label Set 1, Set 2, ..., Set 10.
        subsample_label = "Set " + str(i + 1)


        for molecule_id in MOLECULES_OF_INTEREST:

            print(
                f"{subsample_label:>12} | "
                f"{molecule_id:>12} | "
                f"{subsample_averages[molecule_id][i]:>20.8f} | "
                f"{subsample_uncertainties[molecule_id][i]:>24.8f}"
            )


        # Separate the sub-samples visually.
        print("-" * 82)


############################################################
# Print final molecule analysis results
############################################################

def print_final_results(
    final_results
):

    print()
    print("=" * 75)
    print("FINAL RESULTS - COMPLETE SAMPLE")
    print("=" * 75)


    print(
        f"{'Molecule ID':>12} | "
        f"{'Average per event':>20} | "
        f"{'Statistical uncertainty':>24}"
    )


    print("-" * 75)


    for molecule_id in MOLECULES_OF_INTEREST:

        average, uncertainty = (
            final_results[
                molecule_id
            ]
        )


        print(
            f"{molecule_id:>12} | "
            f"{average:>20.8f} | "
            f"{uncertainty:>24.8f}"
        )


############################################################
# Print overall asymmetry results
############################################################

def print_pair_results(
    pair_results
):

    print()
    print("=" * 90)
    print("OVERALL MOLECULE ASYMMETRY")
    print("=" * 90)


    print(
        f"{'Pair':>20} | "
        f"{'Mean difference':>18} | "
        f"{'Uncertainty':>15} | "
        f"{'|Difference| / uncertainty':>26}"
    )


    print("-" * 90)


    for pair in MOLECULE_PAIRS:

        (
            difference,
            uncertainty,
            significance
        ) = pair_results[pair]


        pair_name = (
            str(pair[0])
            + " vs "
            + str(pair[1])
        )


        print(
            f"{pair_name:>20} | "
            f"{difference:>18.8f} | "
            f"{uncertainty:>15.8f} | "
            f"{significance:>26.2f}"
        )


############################################################
# Print momentum-dependent results
############################################################

def print_momentum_results(
    momentum_results,
    bin_edges
):

    print()
    print("=" * 90)
    print(
        "MOMENTUM-DEPENDENT ASYMMETRY"
    )
    print("=" * 90)


    for pair in MOLECULE_PAIRS:

        (
            differences,
            uncertainties,
            significances
        ) = momentum_results[pair]


        print()
        print(
            "Pair:",
            pair[0],
            "vs",
            pair[1]
        )


        print(
            f"{'Momentum bin':>18} | "
            f"{'Mean difference':>18} | "
            f"{'Uncertainty':>15} | "
            f"{'|Diff.| / uncertainty':>22}"
        )


        print("-" * 82)


        for i in range(
            len(bin_edges) - 1
        ):

            bin_name = (
                f"{bin_edges[i]:.1f} - "
                f"{bin_edges[i + 1]:.1f}"
            )


            print(
                f"{bin_name:>18} | "
                f"{differences[i]:>18.8f} | "
                f"{uncertainties[i]:>15.8f} | "
                f"{significances[i]:>22.2f}"
            )


############################################################
# Plot momentum-dependent asymmetry
############################################################

def plot_momentum_results(
    momentum_results,
    bin_edges
):

    # Centre of every momentum bin.
    bin_centres = (
        bin_edges[:-1]
        + bin_edges[1:]
    ) / 2


    # Half-width of every momentum bin.
    bin_widths = (
        bin_edges[1:]
        - bin_edges[:-1]
    ) / 2


    for pair in MOLECULE_PAIRS:

        (
            differences,
            uncertainties,
            significances
        ) = momentum_results[pair]


        plt.figure(
            figsize=(9, 6)
        )


        # Vertical error bars = statistical uncertainty.
        # Horizontal error bars = momentum-bin width.

        plt.errorbar(
            bin_centres,
            differences,
            xerr=bin_widths,
            yerr=uncertainties,
            fmt="o",
            capsize=4,
            markersize=5
        )


        # Zero means no difference between
        # the two molecule IDs.
        plt.axhline(
            y=0,
            linestyle="--",
            linewidth=1
        )


        plt.xlabel("Momentum P (1e-20 kg m/s)")

        plt.ylabel("Average frequency difference per event")

        plt.title(
            f"Momentum-dependent asymmetry: {pair[0]} vs {pair[1]}"
        )


        plt.grid(
            alpha=0.25
        )


        plt.tight_layout()


        # Save plot for the README.
        plot_filename = (
            "momentum_asymmetry_"
            + str(pair[0])
            + "_"
            + str(pair[1])
            + ".png"
        )


        plt.savefig(
            plot_filename,
            dpi=300
        )


        plt.show()
        plt.close()


############################################################
# Save analysis results to CSV file
############################################################

def save_csv_results(
    subsample_averages,
    subsample_uncertainties,
    events_per_subsample,
    final_results,
    pair_results,
    momentum_results,
    bin_edges
):

    # Create a folder for the CSV files.
    output_folder = "csv_results"
    os.makedirs(output_folder, exist_ok=True)

   ############################################################
    # 1. Individual sub-sample results
    ############################################################

    rows = []

    for i in range(len(events_per_subsample)):

        for molecule_id in MOLECULES_OF_INTEREST:

            rows.append([
                i + 1,
                molecule_id,
                events_per_subsample[i],
                subsample_averages[molecule_id][i],
                subsample_uncertainties[molecule_id][i]
            ])

    with open(
        os.path.join(output_folder, "subsample_results.csv"),
        "w",
        newline=""
    ) as file:

        writer = csv.writer(file)

        writer.writerow([
            "Subsample",
            "Molecule ID",
            "Number of events",
            "Average per event",
            "Poisson uncertainty"
        ])

        writer.writerows(rows)

    ############################################################
    # 2. Final molecule results
    ############################################################

    with open(
        os.path.join(output_folder, "final_molecule_results.csv"),
        "w",
        newline=""
    ) as file:

        writer = csv.writer(file)

        writer.writerow([
            "Molecule ID",
            "Average per event",
            "Subsampling uncertainty"
        ])

        for molecule_id in MOLECULES_OF_INTEREST:

            average, uncertainty = final_results[molecule_id]

            writer.writerow([
                molecule_id,
                average,
                uncertainty
            ])

    ############################################################
    # 3. Overall molecule pair asymmetry
    ############################################################

    with open(
        os.path.join(output_folder, "overall_pair_asymmetry.csv"),
        "w",
        newline=""
    ) as file:

        writer = csv.writer(file)

        writer.writerow([
            "Positive ID",
            "Negative ID",
            "Mean difference",
            "Uncertainty",
            "Difference / uncertainty"
        ])

        for pair in MOLECULE_PAIRS:

            difference, uncertainty, ratio = pair_results[pair]

            writer.writerow([
                pair[0],
                pair[1],
                difference,
                uncertainty,
                ratio
            ])

    ############################################################
    # 4. Momentum results for all pairs
    ############################################################

    momentum_rows = []

    for pair in MOLECULE_PAIRS:

        differences, uncertainties, ratios = momentum_results[pair]

        pair_rows = []

        for i in range(len(bin_edges) - 1):

            row = [
                pair[0],
                pair[1],
                i + 1,
                bin_edges[i],
                bin_edges[i + 1],
                (bin_edges[i] + bin_edges[i + 1]) / 2,
                differences[i],
                uncertainties[i],
                ratios[i]
            ]

            pair_rows.append(row)
            momentum_rows.append(row)

        # Save a separate momentum table for each pair.
        filename = (
            f"momentum_{pair[0]}_{pair[1]}.csv"
        )

        with open(
            os.path.join(output_folder, filename),
            "w",
            newline=""
        ) as file:

            writer = csv.writer(file)

            writer.writerow([
                "Positive ID",
                "Negative ID",
                "Bin number",
                "Momentum minimum",
                "Momentum maximum",
                "Momentum centre",
                "Mean difference",
                "Uncertainty",
                "Difference / uncertainty"
            ])

            writer.writerows(pair_rows)

    # Save all momentum results in one table.
    with open(
        os.path.join(output_folder, "momentum_asymmetry_results.csv"),
        "w",
        newline=""
    ) as file:

        writer = csv.writer(file)

        writer.writerow([
            "Positive ID",
            "Negative ID",
            "Bin number",
            "Momentum minimum",
            "Momentum maximum",
            "Momentum centre",
            "Mean difference",
            "Uncertainty",
            "Difference / uncertainty"
        ])

        writer.writerows(momentum_rows)

    print()
    print("CSV results saved in:", output_folder)


############################################################
# Check that all files exist
############################################################

def check_files(filenames):

    missing_files = []


    for filename in filenames:

        if not os.path.isfile(filename):
            missing_files.append(filename)


    if missing_files:

        print()
        print(
            "The following files could not be found:"
        )


        for filename in missing_files:
            print("-", filename)


        return False


    return True


############################################################
# Main progrmam for all previous functions
############################################################

def main():

    print("=" * 75)
    print("COMPLETE MOLECULE DATA ANALYSIS")
    print("=" * 75)


    print()
    print("Checking data files...")


    # Check files before beginning the long analysis.
    if not check_files(FILENAMES):
        return


    print(
        "All 10 files found."
    )

    print(
        "Starting analysis..."
    )


    try:

        (
            subsample_averages,
            subsample_uncertainties,
            all_momentum_counts,
            events_per_subsample
        ) = analyse_all_subsamples(
            FILENAMES
        )


    except (
        OSError,
        ValueError
    ) as error:

        print()
        print(
            "Analysis stopped:",
            error
        )

        return


    ############################################################
    # QUESTION 1
    # What are the average counts and uncertainties?
    ############################################################

    final_results = calculate_final_results(
        subsample_averages,
        events_per_subsample
    )


    print_subsample_results(
        subsample_averages,
        subsample_uncertainties
    )


    print_final_results(
        final_results
    )


    ############################################################
    # QUESTION 2
    # Is there an overall molecule asymmetry?
    ############################################################

    pair_results = analyse_pairs(
        subsample_averages,
        events_per_subsample
    )


    print_pair_results(
        pair_results
    )


    ############################################################
    # QUESTION 3
    # Is the asymmetry momentum-dependent?
    ############################################################

    momentum_results, bin_edges = analyse_momentum_asymmetry(
        all_momentum_counts,
        events_per_subsample
    )


    print_momentum_results(
        momentum_results,
        bin_edges
    )


    plot_momentum_results(
        momentum_results,
        bin_edges
    )

    # Save the numerical results to CSV files.
    save_csv_results(
        subsample_averages,
        subsample_uncertainties,
        events_per_subsample,
        final_results,
        pair_results,
        momentum_results,
        bin_edges
    )

    print()
    print("=" * 75)

    print(
        "Total number of events analysed:",
        sum(events_per_subsample)
    )

    print("=" * 75)


############################################################
# Run the program
############################################################

if __name__ == "__main__":
    main()

# Python sets __name__ to "__main__" when this file
# is run directly. This makes sure that main() only
# runs when this file is executed, and not if the file
# is imported into another Python program.