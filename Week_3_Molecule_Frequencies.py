##################################################
# Week 3 Deliverables - Molecule Analysis Code 
##################################################
# --> The code reads a molecule data file and selects molecules of interest only 
# --> The code asks user to enter name of a molecule data file themselves 
# --> Then the program gives the complete results table for all molecules of interest

# For each molecule, the program calculates: 
# --> Average number per event and poisson standard deviation = aqrt(average)
# --> Result for Average ± Standard deviation is given


import numpy as np


############################################################
# Define the list of molecules of interest
############################################################

MOLECULES_OF_INTEREST = [
    211, -211,
    321, -321,
    2212, -2212,
    3122, -3122,
    3312, -3312,
    3334, -3334
]


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
# Read one data file
############################################################

def read_file(filename):

    # Initialise count for every molecule ID.
    molecule_counts = {
        molecule_id: 0
        for molecule_id in MOLECULES_OF_INTEREST
    }

    number_of_events = 0
    event_started = False


    with open(filename, "r") as data_file:

        for line in data_file:

            tokens = line.split()

            # Ignore empty lines.
            if not tokens:
                continue


            if is_molecule_line(tokens):

                molecule_id = int(
                    float(tokens[3])
                )

                if molecule_id in molecule_counts:
                    molecule_counts[molecule_id] += 1

                event_started = True


            else:

                # A non-molecule line marks the
                # boundary between events.
                if event_started:
                    number_of_events += 1

                event_started = False


    # Count the final event.
    if event_started:
        number_of_events += 1


    return molecule_counts, number_of_events


############################################################
# Calculate average and uncertainty 
############################################################

def calculate_statistics(total_count, number_of_events):

    # Protect against invalid event counts.
    if number_of_events <= 0:
        raise ValueError("Number of events must be positive.")

    # Average number of molecules per event.
    average = total_count / number_of_events

    # Poisson uncertainty on the average.
    uncertainty = np.sqrt(total_count) / number_of_events

    return average, uncertainty


############################################################
# Print complete table
############################################################

def print_results(
    molecule_counts,
    number_of_events,
    filename
):

    print()
    print("=" * 70)

    # Display the filename so it is clear which
    # sub-sample the results belong to.
    print("MOLECULE RESULTS -", filename)

    print("=" * 70)

    print(
        f"{'Molecule ID':>12} | "
        f"{'Average per event':>20} | "
        f"{'Statistical uncertainty':>24}"
    )

    print("-" * 70)

    for molecule_id in MOLECULES_OF_INTEREST:

        average, uncertainty = calculate_statistics(
            molecule_counts[molecule_id],
            number_of_events
        )

        print(
            f"{molecule_id:>12} | "
            f"{average:>20.8f} | "
            f"{uncertainty:>24.8f}"
        )

    print("=" * 70)


############################################################
# Define the main program for all previous functions
############################################################

def main():

    filename = input(
        "Enter the data filename: "
    ).strip()


    # Protection against empty input.
    if not filename:

        print(
            "Error: no filename was entered."
        )

        return


    print()
    print("Loading file...")


    # Protect against file errors.
    try:

        (
            molecule_counts,
            number_of_events
        ) = read_file(filename)


    except FileNotFoundError:

        print(
            "Error: file not found."
        )

        return


    except PermissionError:

        print(
            "Error: permission denied."
        )

        return


    except OSError as error:

        print(
            "Error while reading file:",
            error
        )

        return


    # Check that events were actually found.
    if number_of_events == 0:

        print(
            "No events were found. "
            "Check the file format."
        )

        return


    print("Loading complete.")

    print(
        "Number of events:",
        number_of_events
    )


    print_results(
        molecule_counts,
        number_of_events,
        filename)

############################################################
# Run the program 
############################################################

if __name__ == "__main__":
    main()

# Python sets __name__ to "__main__" when this file
# is run directly. This makes sure that main() only
# runs when this file is executed, and not if the file
# is imported into another Python program.