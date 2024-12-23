import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def f_exact(x):
    """Standard function f(x) = 1 / (1 + x^2)"""
    return 1.0 / (1.0 + x * x)

def plot_bspline_interpolation(data_dir, output_dir):
    """Plot B-Spline interpolation methods"""

    plt.style.use('seaborn-v0_8')
    colors = ['r--', 'b--', 'g--', 'm--']
    linewidth = 1.5
    legend_fontsize = 14

    x_exact = np.linspace(-5, 5, 1000)
    y_exact = f_exact(x_exact)

    # Plot Cubic B-Spline interpolation
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'k-', label='Actual Function', linewidth=linewidth)

    # Read BSpline3_N11.csv
    file_path = Path(data_dir) / 'BSpline3_N11.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        plt.plot(df['x'], df['BSpline3_NATURAL'], 'r--', label='BSpline3 NATURAL', linewidth=linewidth)
        plt.plot(df['x'], df['BSpline3_COMPLETE'], 'b--', label='BSpline3 COMPLETE', linewidth=linewidth)
        plt.plot(df['x'], df['BSpline3_NOT_A_KNOT'], 'g--', label='BSpline3 NOT_A_KNOT', linewidth=linewidth)
        plt.plot(df['x'], df['BSpline3_SECOND'], 'm--', label='BSpline3 SECOND', linewidth=linewidth)

    plt.title('Cubic B-Spline Interpolation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(fontsize=legend_fontsize)
    plt.grid(True)
    plt.savefig(Path(output_dir) / 'Cubic_BSpline_interpolation.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot Linear B-Spline interpolation
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'k-', label='Actual Function', linewidth=linewidth)

    # Read BSpline1_N11.csv
    file_path = Path(data_dir) / 'BSpline1_N11.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        plt.plot(df['x'], df['BSpline1'], 'b--', label='BSpline1', linewidth=linewidth)

    plt.title('Linear B-Spline Interpolation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(fontsize=legend_fontsize)
    plt.grid(True)
    plt.savefig(Path(output_dir) / 'Linear_BSpline_interpolation.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot Quadratic B-Spline interpolation
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'k-', label='Actual Function', linewidth=linewidth)

    # Read BSpline2_N10.csv
    file_path = Path(data_dir) / 'BSpline2_N10.csv'
    if file_path.exists():
        df = pd.read_csv(file_path)
        plt.plot(df['x'], df['BSpline2'], 'g--', label='BSpline2', linewidth=linewidth)

    plt.title('Quadratic B-Spline Interpolation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(fontsize=legend_fontsize)
    plt.grid(True)
    plt.savefig(Path(output_dir) / 'Quadratic_BSpline_interpolation.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("Starting plot.py...")
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemC'
    output_dir = script_dir.parent.parent / 'figure' / 'problemC'

    data_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_bspline_interpolation(data_dir, output_dir)

if __name__ == "__main__":
    main()
