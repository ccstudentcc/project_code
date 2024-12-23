import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import fsolve

# Define script directory, data directory, and output directory
script_dir = Path(__file__).parent.resolve()
data_dir = script_dir.parent.parent / 'output' / 'intersect'
output_dir = script_dir.parent.parent / 'figure' / 'intersect'

# Print paths for reference
print(f"Data directory: {data_dir}")
print(f"Output directory: {output_dir}")

t = np.linspace(0, 4*np.pi, 200)
a = 0.5
x = 1.5*np.sin(t) + a*np.sin(1.5*t)
y = 1.5*np.cos(t) - a*np.cos(1.5*t)

plt.figure(figsize=(6, 6))
plt.plot(x, y)
plt.title(r"$r_4$ Five-Pointed Star Curve")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# Example correction: Find self-intersection points by solving x(t1)=x(t2), y(t1)=y(t2) with constraints t1<t2, |t1-t2|>threshold
def find_intersections(param_func, t_range, threshold=1e-4):
    intersections = []
    t_values = np.linspace(*t_range, 200)

    def equations(vars):
        t1, t2 = vars
        x1, y1 = param_func(t1)
        x2, y2 = param_func(t2)
        return (x1 - x2, y1 - y2)

    for i in range(len(t_values)):
        for j in range(i+1, len(t_values)):
            guess = (t_values[i], t_values[j])
            sol, infodict, ier, _ = fsolve(equations, guess, full_output=True)
            if ier == 1:
                if (t_range[0] <= sol[0] <= t_range[1] and
                    t_range[0] <= sol[1] <= t_range[1] and
                    abs(sol[0] - sol[1]) > threshold):
                    x_solved, y_solved = param_func(sol[0])
                    # Check for duplicates
                    if not any(np.isclose(x_solved, pt[0], atol=1e-5) and
                               np.isclose(y_solved, pt[1], atol=1e-5) for pt in intersections):
                        intersections.append((x_solved, y_solved))
    return intersections

# Detect self-intersection points using the example curve r4
def parametric_eq(t):
    return 1.5*np.sin(t) + a*np.sin(1.5*t), 1.5*np.cos(t) - a*np.cos(1.5*t)

# Calculate and print self-intersection points
if __name__ == "__main__":
    found_pts = find_intersections(parametric_eq, (0, 4*np.pi))
    print("Intersection points:")
    for pt in found_pts:
        print((round(pt[0], 4), round(pt[1], 4)))

# Plot the intersection points
for pt in found_pts:
    plt.plot(pt[0], pt[1], 'ro')

# Annotate the intersection points
sorted_coords = sorted(found_pts, key=lambda c: (-c[1], c[0]))
for i, coord in enumerate(sorted_coords, 1):
    x_rounded = round(coord[0], 4)
    y_rounded = round(coord[1], 4)
    plt.annotate(f"({x_rounded}, {y_rounded})", (coord[0], coord[1]),
                 textcoords="offset points", xytext=(5,5), ha='left', color='red')

# Save the plot to the output directory
output_file = output_dir / 'five_pointed_star_curve.png'
plt.savefig(output_file)
plt.close()
