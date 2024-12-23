import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 定义保存目录
script_dir = Path(__file__).parent.resolve()
output_dir = script_dir.parent.parent / 'figure' / 'problemA'

# 创建目录（如果不存在）
output_dir.mkdir(parents=True, exist_ok=True)

N = np.array([6, 11, 21, 41, 81])
max_error = np.array([0.421705, 0.020529, 0.003169, 0.000275, 0.000016])

# 计算自然对数
ln_N = np.log(N)
ln_error = np.log(max_error)

# 使用 NumPy 的最小二乘法进行线性拟合
m, b = np.polyfit(ln_N, ln_error, 1)

print(f"Slope: {m}, Intercept: {b}")
print(f"MaxError = exp({b:.2f}) * N^{m:.2f}")

# 绘图
plt.figure(figsize=(8,6))
plt.plot(ln_N, ln_error, 'bo', label='Actual Error')

# 生成拟合线的数据点
ln_N_fit = np.linspace(min(ln_N), max(ln_N), 100)
ln_error_fit = m * ln_N_fit + b
plt.plot(ln_N_fit, ln_error_fit, 'r--', label=f'Least Squares Fit: ln(MaxError)={m:.2f} ln(N)+{b:.2f}')

plt.xlabel('ln(Number of Interpolation Points $N$)')
plt.ylabel('ln(Maximum Error)')

# 保存图形
plt.legend(fontsize=12)
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.tight_layout()
plt.savefig(output_dir / 'lnN_vs_lnMaxError.png', dpi=300)
plt.close()