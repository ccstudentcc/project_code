import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os

def main():
    # Define paths
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemE'
    output_dir = script_dir.parent.parent / 'figure' / 'problemE'
    
    print(f"Script directory: {script_dir}")
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory (if it doesn't exist)
    os.makedirs(output_dir, exist_ok=True)
    
    # Read CSV file
    maxerror_path = data_dir / 'maxerror.csv'
    if not maxerror_path.exists():
        print(f"Error: {maxerror_path} does not exist.")
        return
    
    df = pd.read_csv(maxerror_path)
    
    # Remove data with maxError <= 0 to avoid log errors
    df = df[df['maxError'] > 0]
    
    # Calculate ln(N) and ln(maxError)
    df['ln_N'] = np.log(df['N'])
    df['ln_maxError'] = np.log(df['maxError'])
    
    # Get all curves
    curves = df['curve'].unique()
    
    # Define spline method colors
    spline_colors = {
        'BSpline': 'blue',
        'PPForm': 'red'
    }
    
    # Define boundary condition marker styles
    boundary_markers = {
        'NOT_A_KNOT': 'o',
        'PERIODIC': 's',
        'COMPLETE': '^',
        'SECOND': 'D'
    }
    
    # Iterate through each curve
    for curve in curves:
        plt.figure(figsize=(8, 6))
        curve_data = df[df['curve'] == curve]
        
        spline_types = curve_data['splineType'].unique()
        
        for spline in spline_types:
            spline_data = curve_data[curve_data['splineType'] == spline]
            
            boundary_conditions = spline_data['boundary'].unique()
            
            for boundary in boundary_conditions:
                subset = spline_data[spline_data['boundary'] == boundary]
                
                # Sort data for fitting
                subset = subset.sort_values('N')
                x = subset['ln_N'].values
                y = subset['ln_maxError'].values
                
                # Construct design matrix [x, 1] to compute y = slope * x + intercept
                A = np.vstack([x, np.ones(len(x))]).T
                result = np.linalg.lstsq(A, y, rcond=None)
                slope, intercept = result[0]
                
                # Generate fit line data
                y_fit = slope * x + intercept
                
                # Plot data points
                label = f"{spline}, {boundary}\n$y = {slope:.4f}x + {intercept:.4f}$"
                plt.scatter(x, y, 
                            color=spline_colors.get(spline, 'black'), 
                            marker=boundary_markers.get(boundary, 'o'), 
                            label=label)
                
                # Plot fit line
                plt.plot(x, y_fit, 
                         color=spline_colors.get(spline, 'black'), 
                         linestyle='--')
                
                # Print fit parameters
                # print(f"Curve: {curve}, Spline: {spline}, Boundary: {boundary}")
                # print(f"  Slope: {slope:.4f}, Intercept: {intercept:.4f}")
        
        plt.xlabel('ln(N)', fontsize=12)
        plt.ylabel('ln(maxError)', fontsize=12)
        plt.title(f'Curve: {curve}', fontsize=14)
        plt.legend(fontsize=10, loc='best', markerscale=1.2)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        # Save plot
        output_file = output_dir / f'{curve}_lnMaxError_vs_lnN.png'
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        print(f"Saved plot for curve {curve} at {output_file}\n")
    
    print("All plots have been generated and saved.")

if __name__ == "__main__":
    main()
