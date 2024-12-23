import os
import matplotlib
matplotlib.use('Agg')  # 使用非互动后端
import matplotlib.pyplot as plt
from pathlib import Path
import ast
from mpl_toolkits.mplot3d import Axes3D  # 用于3D绘图

def parse_tuple(s):
    """
    解析字符串形式的元组。
    """
    try:
        return ast.literal_eval(s)
    except Exception as e:
        print(f"Error parsing tuple: {s}\nError: {e}")
        return None

def extract_params(filename):
    """
    从文件名中提取曲线名、N值、Spline类型和节点类型。
    """
    filename_without_ext = os.path.splitext(filename)[0]
    parts = filename_without_ext.split('_')
    
    try:
        N_index = parts.index('N')
        curve_name = '_'.join(parts[:N_index])
        N = parts[N_index + 1]
        splineType = parts[N_index + 2]
        nodeType = parts[N_index + 3]
        return curve_name, N, splineType, nodeType
    except (ValueError, IndexError):
        print(f"Filename {filename} 格式不正确。")
        return None, None, None, None

def main():
    print("Starting plot.py...")
    # 定义输出目录
    script_dir = Path(__file__).parent.resolve()
    data_dir = script_dir.parent.parent / 'output' / 'problemE'
    output_dir = script_dir.parent.parent / 'figure' / 'problemE'

    print(f"Script directory: {script_dir}")
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_dir}")

    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)

    # 获取所有CSV文件，排除maxerror.csv
    csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv') and 'maxerror' not in f]
    print(f"Found {len(csv_files)} CSV files.")

    if not csv_files:
        print("No CSV files found. Exiting.")
        return

    # 创建一个字典来存储所有数据，按N和节点类型分类
    data_dict = {}

    # 记录3D曲线文件
    curve3d_files = []

    for csv_file in csv_files:
        file_path = os.path.join(data_dir, csv_file)
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

        # 跳过表头
        if not lines:
            print(f"{csv_file} 是空文件。")
            continue
        header = lines[0].strip().split(',')
        if header != ['t', 'fittingPoints', 'actualPoints']:
            print(f"{csv_file} 的表头格式不正确。")
            continue

        curve_name, N, splineType, nodeType = extract_params(csv_file)
        if not all([curve_name, N, splineType, nodeType]):
            print(f"Skipping file {csv_file} due to parameter extraction failure.")
            continue  # 跳过无法解析的文件

        print(f"Processing file: {csv_file}")
        print(f"Extracted params - Curve: {curve_name}, N: {N}, SplineType: {splineType}, NodeType: {nodeType}")

        # 初始化data_dict结构
        if N not in data_dict:
            data_dict[N] = {}
        if nodeType not in data_dict[N]:
            data_dict[N][nodeType] = {'r1': [], 'r2': [], 'r3': []}

        for line in lines[1:]:
            parts = []
            current = ''
            in_parentheses = False
            for char in line:
                if char == '(':
                    in_parentheses = True
                    current += char
                elif char == ')':
                    in_parentheses = False
                    current += char
                elif char == ',' and not in_parentheses:
                    parts.append(current.strip())
                    current = ''
                else:
                    current += char
            if current:
                parts.append(current.strip())

            if len(parts) != 3:
                print(f"Skipping line due to incorrect format: {line.strip()}")
                continue

            t_str, fitting_str, actual_str = parts
            t = float(t_str)
            fitting = parse_tuple(fitting_str)
            actual = parse_tuple(actual_str)

            if not isinstance(fitting, tuple) or not isinstance(actual, tuple):
                print(f"Skipping line due to invalid tuple format: {line.strip()}")
                continue

            # 检查是否为3D曲线
            if len(fitting) == 3:
                is_3d = True
            else:
                is_3d = False

            if is_3d:
                curve3d_files.append(csv_file)
                print(f"Identified {csv_file} as 3D curve.")
                break  # 之后单独处理3D曲线
            else:
                # 根据curve_name分类
                curve_key = None
                if curve_name.lower().startswith('r1'):
                    curve_key = 'r1'
                elif curve_name.lower().startswith('r2'):
                    curve_key = 'r2'
                elif curve_name.lower().startswith('r3'):
                    curve_key = 'r3'
                else:
                    print(f"Unknown curve name {curve_name} in file {csv_file}. Skipping.")
                    break

                entry = {
                    'fittingPoints': fitting,
                    'actualPoints': actual
                }
                data_dict[N][nodeType][curve_key].append(entry)

        if is_3d:
            curve3d_files.append(csv_file)

    print("Done processing CSV files.")
    print(f"Preparing to plot data for {len(data_dict)} different N values.")

    # 打印 data_dict 结构以确认数据分类情况
    print("Final data_dict structure:")
    for N, node_types in data_dict.items():
        for nodeType, curves in node_types.items():
            print(f"N={N}, NodeType={nodeType}: r1={len(curves['r1'])}, r2={len(curves['r2'])}, r3={len(curves['r3'])}")

    # 开始绘图
    for N in sorted(data_dict.keys(), key=lambda x: int(x)):
        for nodeType in data_dict[N]:
            print(f"Plotting for N={N}, NodeType={nodeType}...")
            plt.figure(figsize=(10, 6))
            plt.title(f'N={N}, NodeType={nodeType}')
            plt.xlabel('X')
            plt.ylabel('Y')

            # 记录是否已经绘制过原曲线
            original_plot_done = False

            # 绘制r1, r2, r3的拟合曲线
            for curve in ['r1', 'r2', 'r3']:
                entries = data_dict[N][nodeType][curve]
                if not entries:
                    print(f"No entries for {curve.upper()} in N={N}, NodeType={nodeType}.")
                    continue
                for entry in entries:
                    fitting = entry['fittingPoints']
                    actual = entry['actualPoints']

                    # 确保 fitting 和 actual 是元组
                    if not isinstance(fitting, tuple) or not isinstance(actual, tuple):
                        print(f"Invalid data types for N={N}, NodeType={nodeType}, Curve={curve.upper()}")
                        continue

                    # 分离坐标
                    try:
                        actual_x = [actual[0]]
                        actual_y = [actual[1]]
                        fitting_x = [fitting[0]]
                        fitting_y = [fitting[1]]
                    except TypeError as te:
                        print(f"TypeError while processing points for N={N}, NodeType={nodeType}, Curve={curve.upper()}: {te}")
                        continue

                    if not original_plot_done:
                        plt.plot(actual_x, actual_y, 'k-', label='Original Curve')
                        original_plot_done = True

                    label = f"{curve.upper()} | {splineType} | {nodeType}"
                    plt.plot(fitting_x, fitting_y, 'b--', label=label)

            plt.legend(loc='best', fontsize='small')
            plt.grid(True)

            # 保存图形到文件
            output_file_path = os.path.join(output_dir, f"N_{N}_{nodeType}.png")
            plt.savefig(output_file_path, dpi=300)
            plt.close()
            print(f"Saved 2D plot: {output_file_path}")

    print("Finished plotting 2D curves.")

    # 处理3D曲线
    for csv_file in curve3d_files:
        print(f"Plotting 3D curve for file: {csv_file}")
        file_path = os.path.join(data_dir, csv_file)
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue

        # 跳过表头
        if not lines:
            print(f"{csv_file} 是空文件。")
            continue
        header = lines[0].strip().split(',')
        if header != ['t', 'fittingPoints', 'actualPoints']:
            print(f"{csv_file} 的表头格式不正确。")
            continue

        curve_name, N, splineType, nodeType = extract_params(csv_file)
        if not all([curve_name, N, splineType, nodeType]):
            print(f"Skipping file {csv_file} due to parameter extraction failure.")
            continue  # 跳过无法解析的文件

        print(f"Processing 3D file: {csv_file}")
        
        # 初始化data_dict结构
        if N not in data_dict:
            data_dict[N] = {}
        if nodeType not in data_dict[N]:
            data_dict[N][nodeType] = {'r1': [], 'r2': [], 'r3': []}

        for line in lines[1:]:
            parts = []
            current = ''
            in_parentheses = False
            for char in line:
                if char == '(':
                    in_parentheses = True
                    current += char
                elif char == ')':
                    in_parentheses = False
                    current += char
                elif char == ',' and not in_parentheses:
                    parts.append(current.strip())
                    current = ''
                else:
                    current += char
            if current:
                parts.append(current.strip())

            if len(parts) != 3:
                print(f"Skipping line due to incorrect format: {line.strip()}")
                continue

            t_str, fitting_str, actual_str = parts
            t = float(t_str)
            fitting = parse_tuple(fitting_str)
            actual = parse_tuple(actual_str)

            if not isinstance(fitting, tuple) or not isinstance(actual, tuple):
                print(f"Skipping line due to invalid tuple format: {line.strip()}")
                continue

            entry = {
                'fittingPoints': fitting,
                'actualPoints': actual
            }

            # 根据curve_name分类
            curve_key = None
            if curve_name.lower().startswith('r1'):
                curve_key = 'r1'
            elif curve_name.lower().startswith('r2'):
                curve_key = 'r2'
            elif curve_name.lower().startswith('r3'):
                curve_key = 'r3'
            else:
                print(f"Unknown curve name {curve_name} in file {csv_file}. Skipping.")
                continue

            data_dict[N][nodeType][curve_key].append(entry)

        # 绘制3D曲线
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f'N={N}, NodeType={nodeType}')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        original_plot_done_3d = False

        for curve in ['r1', 'r2', 'r3']:
            entries = data_dict[N][nodeType][curve]
            if not entries:
                print(f"No entries for {curve.upper()} in N={N}, NodeType={nodeType}.")
                continue
            for entry in entries:
                fitting = entry['fittingPoints']
                actual = entry['actualPoints']

                # 确保 fitting 和 actual 是元组
                if not isinstance(fitting, tuple) or not isinstance(actual, tuple):
                    print(f"Invalid data types for N={N}, NodeType={nodeType}, Curve={curve.upper()}")
                    continue

                # 分离坐标
                try:
                    actual_x, actual_y, actual_z = actual
                    fitting_x, fitting_y, fitting_z = fitting
                except TypeError as te:
                    print(f"TypeError while processing points for N={N}, NodeType={nodeType}, Curve={curve.upper()}: {te}")
                    continue

                if not original_plot_done_3d:
                    ax.plot(actual_x, actual_y, actual_z, 'k-', label='Original Curve')
                    original_plot_done_3d = True

                label = f"{curve.upper()} | {splineType} | {nodeType}"
                ax.plot(fitting_x, fitting_y, fitting_z, 'b--', label=label)

        ax.legend(loc='best', fontsize='small')
        ax.grid(True)

        # 保存图形到文件
        output_file_path = os.path.join(output_dir, f"N_{N}_{nodeType}_3D.png")
        plt.savefig(output_file_path, dpi=300)
        plt.close()
        print(f"Saved 3D plot: {output_file_path}")

    print("Finished plotting 3D curves.")
    print("All plots have been generated and saved successfully.")

if __name__ == "__main__":
    main()