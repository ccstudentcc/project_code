import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def plot_bspline_curves():
    # Define script directory, data directory, and output directory
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'requirement5'
    output_dir = script_dir.parent.parent / 'figure' / 'requirement5'
    
    # Print paths for reference
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory (if it doesn't exist)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define files to read
    bspline_files = {
        'Cubic B-spline': 'cubic_bspline.csv',
        'Quintic B-spline': 'quintic_bspline.csv'
    }
    
    plt.figure(figsize=(10, 6))
    
    # Colors and styles list
    colors = ['r-', 'b--', 'g-.', 'm:']
    
    for idx, (label, filename) in enumerate(bspline_files.items()):
        file_path = data_dir / filename
        if file_path.exists():
            # Read CSV file
            data = pd.read_csv(file_path, header=None, names=['x', 'y'])
            x = data['x']
            y = data['y']
            
            # Plot curve
            plt.plot(x, y, colors[idx % len(colors)], label=label)
            print(f"Plotted file: {file_path}")
        else:
            print(f"File not found: {file_path}")
    
    # Plot original function (if needed)
    # For example, plot f_exact(x) = 1 / (1 + 25x^2)
    # x_exact = np.linspace(-1, 1, 1000)
    # y_exact = 1.0 / (1.0 + 25.0 * x_exact**2)
    # plt.plot(x_exact, y_exact, 'k-', label='Exact Function')
    
    # Set plot title and labels
    plt.title('B-spline Curve Plotting')
    plt.xlabel('x')
    plt.ylabel('y')
    
    # Add legend
    plt.legend(fontsize=12)
    
    # Add grid
    plt.grid(True, linestyle='--', linewidth=0.5)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    output_file = output_dir / 'bspline_curves.png'
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to: {output_file}")
    
    # # Show plot
    # plt.show()
    plt.close()

if __name__ == "__main__":
    plot_bspline_curves()