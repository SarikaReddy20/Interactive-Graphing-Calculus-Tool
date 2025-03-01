import tkinter as tk
from tkinter import simpledialog,messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff, integrate, lambdify, solve, sin, cos, tan, cot, sec, csc, pi, log, exp, Interval

class MathPlotApp:

    def __init__(self, master):
        self.master = master
        master.title("MathPlot")
        self.option_label = tk.Label(master, text="Choose an option:")
        self.option_label.pack()
        self.intersection_button = tk.Button(master, text="Intersection Graph", command=self.plot_intersection)
        self.intersection_button.pack()
        self.differentiation_button = tk.Button(master, text="Differentiation Graph", command=self.plot_differentiation)
        self.differentiation_button.pack()
        self.integration_button = tk.Button(master, text="Integration Graph", command=self.plot_integration)
        self.integration_button.pack()
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.canvas.get_tk_widget().pack()

    def plot_intersection(self):
        self.ax.clear()
        equations = []
        num_equations = simpledialog.askinteger("Input", "Enter the number of equations:")
        for i in range(num_equations):
            equation = simpledialog.askstring("Input", f"Enter equation {i + 1}:")
            equations.append(equation)
            # Find and display roots of the equation
            roots = self.find_roots(equation)
            if roots:
                print(f"Roots of Equation {i + 1}: {roots}")
            else:
                print(f"No real roots found for Equation {i + 1}")

        x_min = simpledialog.askfloat("Input", "Enter the minimum value of x:")
        x_max = simpledialog.askfloat("Input", "Enter the maximum value of x:")
        self.plot_intersection_graph(equations, x_min, x_max)

    def plot_intersection_graph(self, equations, x_min, x_max):
        values_of_x = np.arange(x_min, x_max + 0.1, 0.1)
        intersection_points = self.find_intersection_points(equations, values_of_x)
        for i, equation in enumerate(equations):
            values_of_y = [self._y_function(x, equation) for x in values_of_x]
            self.ax.plot(values_of_x, values_of_y, label=f'Equation {i + 1}')
        # Annotate and plot intersection points
        for point in intersection_points:
            self.ax.scatter(point[0], point[1], color='red', marker='o')
            self.ax.annotate(f'({point[0]:.2f}, {point[1]:.2f})', (point[0], point[1]),
                             textcoords="offset points", xytext=(5, -5), ha='center')
        self.ax.axhline(0, color='black', linewidth=0.5)  # Horizontal axis y=0
        self.ax.axvline(0, color='black', linewidth=0.5)  # Vertical axis x=0
        self.ax.legend()
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.grid(True)  # Show grid lines
        self.canvas.draw()

    def plot_differentiation(self):
        self.ax.clear()
        equation = simpledialog.askstring("Input", 'Enter the equation to differentiate (e.g., "x**3 - 8*x + 9"):')
        start = simpledialog.askfloat("Input", 'Enter the start of the range:')
        end = simpledialog.askfloat("Input", 'Enter the end of the range:')
        color = simpledialog.askstring("Input", 'Enter the color of the graph (e.g., "red", "blue"):')
        # Find and display roots of the equation
        roots = self.find_roots(equation)
        if roots:
            print(f"Roots of the Equation: {roots}")
        else:
            print("No real roots found for the Equation")
        self.plot_differentiation_graph(equation, start, end, color)

    def plot_differentiation_graph(self, equation, start, end, color):
        # Define the variable
        x = symbols('x')
        # Parse the equation
        expr = eval(equation)
        # Compute the first and second derivatives
        expr_diff1 = diff(expr, x)
        expr_diff2 = diff(expr_diff1, x)
        # Convert expressions to lambda functions
        f = lambdify(x, expr, modules=['numpy'])
        f_diff1 = lambdify(x, expr_diff1, modules=['numpy'])
        f_diff2 = lambdify(x, expr_diff2, modules=['numpy'])
        # Generate x values for plotting
        x_values = np.linspace(start, end, 1000)
        # Plot the equation and its derivatives
        self.ax.plot(x_values, f(x_values), label='Equation', color=color)
        self.ax.plot(x_values, f_diff1(x_values), label="1st Derivative", linestyle='--', color=color)
        self.ax.plot(x_values, f_diff2(x_values), label="2nd Derivative", linestyle=':', color=color)
        # Add labels and legend
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_title('Graph of Equation and Its Derivatives')
        self.ax.legend()
        self.ax.grid(True)
        # Show the plot
        self.canvas.draw()

    def plot_integration(self):
        self.ax.clear()
        equation = simpledialog.askstring("Input", 'Enter the equation to integrate (e.g., "x**3 - 8*x + 9"):')
        start = simpledialog.askfloat("Input", 'Enter the start of the range:')
        end = simpledialog.askfloat("Input", 'Enter the end of the range:')
        color = simpledialog.askstring("Input", 'Enter the color of the graph (e.g., "red", "blue"):')
        # Find and display roots of the equation
        roots = self.find_roots(equation)
        if roots:
            print(f"Roots of the Equation: {roots}")
        else:
            print("No real roots found for the Equation")
        self.plot_integration_graph(equation, start, end, color)

    def plot_integration_graph(self, equation, start, end, color):
        # Define the variable
        x = symbols('x')
        # Parse the equation
        expr = eval(equation)
        # Integrate the expression
        integral_expr = integrate(expr, x)
        # Convert expressions to lambda functions
        f = lambdify(x, expr, modules=['numpy'])
        F = lambdify(x, integral_expr, modules=['numpy'])
        # Generate x values for plotting
        x_values = np.linspace(start, end, 1000)
        # Plot the equation and its integral
        self.ax.plot(x_values, f(x_values), label='Equation', color=color)
        self.ax.plot(x_values, F(x_values), label="Integral", linestyle='--', color=color)
        # Add labels and legend
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_title('Graph of Equation and Its Integral')
        self.ax.legend()
        self.ax.grid(True)
        # Show the plot
        self.canvas.draw()
    def _y_function(self, x, equation):
        try:
            return eval(equation)
        except Exception:
            return np.nan  # Return NaN if there's an error in the evaluation

    def find_intersection_points(self, equations, x_values):
        intersection_points = []
        for equation1 in equations:
            for equation2 in equations:
                if equation1 != equation2:  # Avoid comparing the same equation
                    y_values1 = np.array([self._y_function(x, equation1) for x in x_values])
                    y_values2 = np.array([self._y_function(x, equation2) for x in x_values])
                    if y_values1.dtype != 'float' or y_values2.dtype != 'float':
                        continue  # Skip non-numeric values
                    intersection_indices = np.where(np.isclose(y_values1, y_values2, atol=0.01))[0]
                    for index in intersection_indices:
                        intersection_points.append((x_values[index], y_values1[index]))
        return intersection_points

    def find_roots(self, equation):
        x = symbols('x')
        expr = eval(equation, {"x": x, "sin": sin, "cos": cos, "tan": tan, "cot": cot, "sec": sec, "csc": csc, "pi": pi,
                           "log": log, "exp": exp})
        try:
            roots = solve(expr, x)
            real_roots = [root.evalf() for root in roots if root.is_real]
            return real_roots
        except NotImplementedError:
            # For equations with multiple generators, try to find numerical approximations using a numerical method
            # Here, we use the bisection method for simplicity
            x_min, x_max = 0.1, 10  # Define a range for root searching
            roots = []
            while x_min < x_max:
                try:
                    root = solve(expr, x, domain=Interval(x_min, x_max))
                    if root:
                        roots.append(root[0].evalf())
                    x_min = x_max
                except:
                    pass
                x_min += 0.1
            return roots
root = tk.Tk()
app = MathPlotApp(root)
root.mainloop()