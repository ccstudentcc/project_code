# plotA.py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def f_exact(x):
    """Exact function f(x) = 1/(1+25x^2)"""
    return 1.0 / (1.0 + 25.0 * x * x)

def plot_interpolation_methods(N, output_dir, data_dir):
    """Plot interpolation results for a given N"""
    
    plt.style.use('seaborn-v0_8')
    colors = ['r--', 'b--', 'g--', 'm--']
    linewidth = 1.5
    legend_fontsize = 14  
    
    # Generate exact function curve
    x_exact = np.linspace(-1, 1, 1000)
    y_exact = f_exact(x_exact)
    
    # 1. Plot PPForm3
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'k-', label='Exact', linewidth=linewidth)
    
    # Read PPForm3 data
    file_path = Path(data_dir) / f'PPForm3_N{N}.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        boundary_conditions = ['NATURAL', 'COMPLETE', 'NOT_A_KNOT', 'SECOND']
        for bc, color in zip(boundary_conditions, colors):
            plt.plot(df['x'], df[f'PPForm3_{bc}'], color, 
                    label=f'PPForm3 {bc}', linewidth=linewidth)
    
    plt.title(f'PPForm Cubic Spline Interpolation (N={N})')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(fontsize=legend_fontsize)
    plt.grid(True)
    plt.savefig(Path(output_dir) / f'PPForm3_N{N}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Plot BSpline3
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'k-', label='Exact', linewidth=linewidth)
    
    # Read BSpline3 data
    file_path = Path(data_dir) / f'BSpline3_N{N}.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        boundary_conditions = ['NATURAL', 'COMPLETE', 'NOT_A_KNOT', 'SECOND']
        for bc, color in zip(boundary_conditions, colors):
            plt.plot(df['x'], df[f'BSpline3_{bc}'], color, 
                    label=f'BSpline3 {bc}', linewidth=linewidth)
    
    plt.title(f'BSpline Cubic Spline Interpolation (N={N})')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(fontsize=legend_fontsize)
    plt.grid(True)
    plt.savefig(Path(output_dir) / f'BSpline3_N{N}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Plot linear interpolation methods
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'k-', label='Exact', linewidth=linewidth)
    
    # Plot PPForm1
    file_path = Path(data_dir) / f'PPForm1_N{N}.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        plt.plot(df['x'], df['PPForm1'], 'r--', 
                label='PPForm1', linewidth=linewidth)
    
    # Plot BSpline1
    file_path = Path(data_dir) / f'BSpline1_N{N}.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        plt.plot(df['x'], df['BSpline1'], 'b--', 
                label='BSpline1', linewidth=linewidth)
    
    plt.title(f'Linear Spline Interpolation (N={N})')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(fontsize=legend_fontsize+2)
    plt.grid(True)
    plt.savefig(Path(output_dir) / f'linear_N{N}.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemA'
    output_dir = script_dir.parent.parent / 'figure' / 'problemA'
    
    data_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    N_values = [6, 11, 21, 41, 81]
    for N in N_values:
        print(f"Processing plots for N = {N}")
        plot_interpolation_methods(N, output_dir, data_dir)

if __name__ == "__main__":
    main()