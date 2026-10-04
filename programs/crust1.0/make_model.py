import sys
import crust1

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: make_model.py evlo evla")
    sys.exit(1)

# Parse the command-line arguments
lat = float(sys.argv[2])
lon = float(sys.argv[1])

# Make model instance
model = crust1.crustModel()

# Run the model for each lat and lon pair
model_result = model.get_point(lat, lon)

# Save the velocity model to an ASCII file
with open('vmodel.txt', 'w') as f:
    for i, layer in enumerate(model_result):
        values = model_result[layer]
        
        # Round values to two decimal places and exclude the last column
        values = [round(v, 2) for v in values[:-1]]
        
        # Add 0.01 to the second and fourth columns of the first row
        if i == 0:
            values[1] += 0.01
            values[3] += 0.01
        
        # Write the processed values to the file
        f.write(' '.join(f"{v:.2f}" for v in values) + '\n')
print("1D velocity model saved to vmodel.txt")
