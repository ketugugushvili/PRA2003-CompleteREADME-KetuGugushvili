#######################################################################
# Week 2: Momentum 
#######################################################################

# Start with making a header for the output table
# # Display id and momentum as the header labels for the 2 columns 
print(f"\n{'ID':<6} | {'Momentum':>10}")
print("-" * 20)

# Open the input file for the program to read the molecule data 
with open("output-Set0.txt", "r") as infile:
    # # The code should be able to open "output-Set0.txt" file from any computer if they are in the same folder
    # # # If the input file still doesn't open, copy the file path and paste it inside "" instead of the file name
    
    # Skip the event header (first line) of the data file 
    next(infile) 

    # Read remaining data lines for each molecule 
    for line in infile:

        # Split the line into its four values
        # # three momentum components and the molecule id
        px, py, pz, molecule_id = line.split()

        # Convert the values from text to numbers 
        px = float(px)
        py = float(py)
        pz = float(pz)
        molecule_id = int(molecule_id)

        # Import math to use built-in math functions 
        import math 
        # Define a function to calculate total momentum for each molecule with its formula
        # # Use x, y, and z momentum components 
        def totalmomentum(px, py, pz):
            momentum = math.sqrt(px**2 + py**2 + pz**2)
            return momentum 

        # Call the function that was defined above to calculate the total momentum
        momentum = totalmomentum(px, py, pz)

        # Print id and momentum values under the header for each molecule
        print(f"{molecule_id:<6} | {momentum:>10.4f}")