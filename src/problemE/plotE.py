import os
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import ast
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting
from multiprocessing import Pool, cpu_count

def parse_point(line):
    try:
        parts = line.strip().split(',')
        if len(parts) not in [2, 3]:
            print(f"Invalid number of coordinates in line: {line}")
            return None
        point = tuple(float(part) for part in parts)
        return point
    except Exception as e:
        print(f"Error parsing line: {line}\nError: {e}")
        return None

def extract_params(filename):
    filename_without_ext = os.path.splitext(filename)[0]
    parts = filename_without_ext.split('_')
    
    if parts[-1] == 'fitted':
        try:
            curve_name = parts[0]
            N = parts[2]
            splineType = parts[3]
            nodeType = parts[4]
            boundary = parts[5]
            return curve_name, N, splineType, nodeType, boundary
        except IndexError:
            print(f"Filename {filename} is missing parameters.")
            return None
    elif parts[-1] == 'actual':
        try:
            curve_name = parts[0]
            return curve_name, None, None, None, None
        except IndexError:
            print(f"Filename {filename} is missing curve name.")
            return None
    else:
        print(f"Unknown file type for filename: {filename}")
        return None

def load_points(file_path):
    points = []
    with open(file_path, 'r') as f:
        for line in f:
            point = parse_point(line)
            if point:
                points.append(point)
    return points

def plot_2d(fitted_points, actual_points, curve_name, N, splineType, nodeType, boundary, output_dir):
    plt.figure(figsize=(10, 6))
    plt.title(f'Curve: {curve_name}, N={N}, Spline={splineType}, NodeType={nodeType}, Boundary={boundary}', fontsize=16)
    plt.xlabel('X', fontsize=14)
    plt.ylabel('Y', fontsize=14)

    actual_x = [point[0] for point in actual_points]
    actual_y = [point[1] for point in actual_points]
    plt.plot(actual_x, actual_y, 'k-', label='Actual Curve')

    fitted_x = [point[0] for point in fitted_points]
    fitted_y = [point[1] for point in fitted_points]
    plt.plot(fitted_x, fitted_y, 'b--', label='Fitted Curve')

    plt.legend(loc='best', fontsize='medium')
    plt.grid(True)

    output_file_path = os.path.join(output_dir, f"N_{N}_{curve_name}_{splineType}_{nodeType}_{boundary}.png")
    plt.savefig(output_file_path, dpi=300)
    plt.close()

def plot_3d(fitted_points, actual_points, curve_name, N, splineType, nodeType, boundary, output_dir):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Curve: {curve_name}, N={N}, Spline={splineType}, NodeType={nodeType}, Boundary={boundary}', fontsize=16)
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_zlabel('Z', fontsize=14)

    actual_x = [point[0] for point in actual_points]
    actual_y = [point[1] for point in actual_points]
    actual_z = [point[2] for point in actual_points]
    ax.plot(actual_x, actual_y, actual_z, 'k-', label='Actual Curve')

    fitted_x = [point[0] for point in fitted_points]
    fitted_y = [point[1] for point in fitted_points]
    fitted_z = [point[2] for point in fitted_points]
    ax.plot(fitted_x, fitted_y, fitted_z, 'b--', label='Fitted Curve')

    ax.legend(loc='best', fontsize='medium')
    ax.grid(True)

    output_file_path = os.path.join(output_dir, f"N_{N}_{curve_name}_{splineType}_{nodeType}_{boundary}_3D.png")
    plt.savefig(output_file_path, dpi=300)
    plt.close()

def process_file(fitted_file, actual_dict, data_dir, output_dir):
    params = extract_params(fitted_file)
    if not params:
        print(f"Skipping file {fitted_file} due to parameter extraction failure.")
        return
    curve_name, N, splineType, nodeType, boundary = params

    actual_points = actual_dict.get(curve_name)
    if not actual_points:
        print(f"No actual points found for curve: {curve_name}. Skipping {fitted_file}.")
        return

    fitted_path = os.path.join(data_dir, fitted_file)
    fitted_points = load_points(fitted_path)
    if not fitted_points:
        print(f"No fitted points found in {fitted_file}. Skipping.")
        return

    is_3d = len(fitted_points[0]) == 3

    if is_3d:
        plot_3d(fitted_points, actual_points, curve_name, N, splineType, nodeType, boundary, output_dir)
    else:
        plot_2d(fitted_points, actual_points, curve_name, N, splineType, nodeType, boundary, output_dir)

def main():
    print("Starting plot.py...")
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemE'
    output_dir = script_dir.parent.parent / 'figure' / 'problemE'

    print(f"Script directory: {script_dir}")
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")

    os.makedirs(output_dir, exist_ok=True)

    all_files = os.listdir(data_dir)
    fitted_files = [f for f in all_files if f.endswith('_fitted.csv')]
    actual_files = [f for f in all_files if f.endswith('_actual.csv')]

    print(f"Found {len(fitted_files)} fitted CSV files and {len(actual_files)} actual CSV files.")

    actual_dict = {}
    for actual_file in actual_files:
        params = extract_params(actual_file)
        if params:
            curve_name = params[0]
            actual_path = os.path.join(data_dir, actual_file)
            actual_points = load_points(actual_path)
            actual_dict[curve_name] = actual_points
            print(f"Loaded actual points for curve: {curve_name} from {actual_file}")

    with Pool(cpu_count()) as pool:
        pool.starmap(process_file, [(fitted_file, actual_dict, data_dir, output_dir) for fitted_file in fitted_files])

    print("All plots have been generated and saved successfully.")

if __name__ == "__main__":
    main()
