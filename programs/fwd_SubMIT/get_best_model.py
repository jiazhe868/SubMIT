import glob

def find_best_model():
    # Initialize variables to store the best model information
    min_misfit = float('inf')
    best_model = None

    # Take ONLY the highest-priority mode file present. Never min-compare
    # across mode files: a stale best_model_exploration.dat from an OLD run
    # (different objective era) carries an incomparable misfit and silently
    # wins, so the fwd renders a month-old model (Calama 3sub, 2026-08-14).
    import os
    best_files = []
    for f in ('best_model_hybrid.dat', 'best_model_ensemble.dat',
              'best_model_exploration.dat'):
        if os.path.exists(f):
            best_files = [f]
            break
    if not best_files:
        best_files = glob.glob('*best.dat')

    # Iterate through each file
    for best_file in best_files:
        with open(best_file, 'r') as f:
            # Read the first line of the file
            line = f.readline().strip().split()

            # Parse the line to get nsub, misfit, and source parameters
            nsub = int(line[0])  # number of subevents
            misfit = float(line[2])  # misfit value

            # Check if this model has a lower misfit than the current best
            if misfit < min_misfit:
                min_misfit = misfit
                best_model = line[2:]  # Store the source parameters (skip first 3 columns)

    # Output the best model's source parameters to Input.model
    if best_model:
        with open('Input.model', 'w') as f:
            for i in range(len(best_model) // 8):  # stride 8: misfit + 7 params per subevent
                # Write each subevent's parameters on a new line
                subevent_params = best_model[i*8+1:(i+1)*8]
                f.write(' '.join(subevent_params) + '\n')
        with open('misfit.dat', 'w') as f:
            f.write(f"{min_misfit}\n")

        print(f"Best model found with misfit {min_misfit} and saved to Input.model")
    else:
        print("No best model files found or files are empty.")

if __name__ == "__main__":
    find_best_model()

