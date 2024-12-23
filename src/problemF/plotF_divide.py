import sys
import re
import matplotlib.pyplot as plt
from pathlib import Path

def parse_divided_truncate_file(file_path):
    """
    Parse the divided_X.csv file to extract (i, j) and their corresponding (x, y) data points.
    
    Parameters:
        file_path (Path): Path to the CSV file.
    
    Returns:
        dict: Keys are tuples (i, j), values are lists of (x, y) data points.
    """
    divided_data = {}
    current_key = None
    pattern = re.compile(r'\((\d+),\s*(\d+)\)')
    
    with file_path.open('r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            match = pattern.match(line)
            if match:
                i, j = int(match.group(1)), int(match.group(2))
                current_key = (i, j)
                divided_data[current_key] = []
                # print(f"Parsed new plot: (i={i}, j={j})")
            else:
                if current_key is not None:
                    try:
                        x, y = map(float, line.split(','))
                        divided_data[current_key].append((x, y))
                    except ValueError:
                        print(f"Warning: Unable to parse line: {line}")
                else:
                    print(f"Warning: Data point found without specified (i,j): {line}")
    return divided_data

def plot_divided_truncate(divided_data, output_file):
    """
    Generate plots based on the parsed divided_data.
    
    Parameters:
        divided_data (dict): Keys are tuples (i, j), values are lists of (x, y) data points.
        output_file (Path): Path to the output image file.
    """
    if not divided_data:
        print("No data to plot.")
        return

    # Determine the maximum i and j values
    max_i = max(key[0] for key in divided_data.keys())
    max_j = max(key[1] for key in divided_data.keys())

    # print(f"Maximum i value: {max_i}, Maximum j value: {max_j}")

    # Create a grid of subplots
    fig, axes = plt.subplots(max_i + 1, max_j + 1, figsize=(4 * (max_j + 1), 3 * (max_i + 1)),
                             squeeze=False)

    # Iterate over all possible (i, j) and plot
    for i in range(max_i + 1):
        for j in range(max_j + 1):
            ax = axes[i][j]
            if j <= i:
                key = (i, j)
                if key in divided_data:
                    x, y = zip(*divided_data[key])
                    ax.plot(x, y, 'b-')  # Blue solid line
                    ax.set_title(f'({i}, {j})')
                else:
                    print(f"Warning: No data found for (i={i}, j={j})")
                # Remove y-axis labels
                ax.set_yticks([])
            else:
                # Hide unnecessary subplots
                ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Figure saved to: {output_file}")
    plt.close()

def plot_divide_truncate_from_file(csv_path):
    # Define script directory, data directory, and output directory
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemF'
    output_dir = script_dir.parent.parent / 'figure' / 'problemF'
    
    # Print paths for reference
    # print(f"Data directory: {data_dir}")
    # print(f"Output directory: {output_dir}")
    
    # Create output directory (if it doesn't exist)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    input_file = Path(csv_path).resolve()
    
    if not input_file.exists():
        print(f"Error: File does not exist - {input_file}")
        sys.exit(1)
    
    # Parse the CSV file
    divided_data = parse_divided_truncate_file(input_file)
    
    # Define output file name
    input_stem = input_file.stem  # e.g., 'divided_2'
    output_file = output_dir / f'{input_stem}.png'
    
    # Plot and save the figure
    plot_divided_truncate(divided_data, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plotF_divide.py <path_to_divided_X.csv>")
        sys.exit(1)
    
    csv_path = sys.argv[1]
    plot_divide_truncate_from_file(csv_path)
