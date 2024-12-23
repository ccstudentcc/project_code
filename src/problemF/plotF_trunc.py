import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def plot_truncate(csv_filename):
    # Define script directory, data directory, and output directory
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemF'
    output_dir = script_dir.parent.parent / 'figure' / 'problemF'
    
    # Print paths for reference
    # print(f"Data directory: {data_dir}")
    # print(f"Output directory: {output_dir}")
    
    # Create output directory (if it doesn't exist)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define input file path
    input_file = data_dir / csv_filename
    
    if not input_file.exists():
        print(f"Error: File does not exist - {input_file}")
        sys.exit(1)
    
    # Read CSV file
    try:
        data = pd.read_csv(input_file, header=None, names=['x', 'y'])
        x = data['x']
        y = data['y']
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    # Plot the graph
    plt.figure(figsize=(8,6))
    plt.plot(x, y, 'b-', label='Truncate Function')
    
    # Set graph title and labels
    plt.title('Truncate Function Plot')
    plt.xlabel('x')
    plt.ylabel('y')
    
    # Add legend
    plt.legend(fontsize=12)
    
    # Add grid
    plt.grid(True, linestyle='--', linewidth=0.5)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the graph
    output_file = output_dir / f'{input_file.stem}.png'
    plt.savefig(output_file, dpi=300)
    print(f"Figure saved to: {output_file}")
    
    # Close the plot
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plotF_trunc.py <csv_filename>")
        sys.exit(1)
    
    csv_filename = sys.argv[1]
    plot_truncate(csv_filename)