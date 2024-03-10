import matplotlib.pyplot as plt
import csv
import numpy as np

def plot_points_from_file(filename):
    x_points, y_points = [], []

    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            x_points.append(float(row[0]))
            y_points.append(float(row[1]))

    fig, ax = plt.subplots()
    ax.scatter(x_points, y_points)

    # Plot unit circles around each point
    for x, y in zip(x_points, y_points):
        circle = plt.Circle((x, y), 1, color='red', fill=False)  # Adjust the radius (10 in this case)
        ax.add_patch(circle)

    ax.set_aspect('equal', adjustable='box')
    plt.title('n = 55, r = 1.045816141550128')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

# Example usage
plot_points_from_file('Centers_before.csv')
plot_points_from_file('Centers_after.csv')
