import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from scipy.optimize import fsolve


def lissage(signal_brut, L):
    # Pad the signal at both ends with reflected values
    padded_signal = np.pad(signal_brut, (L, L), mode='reflect')
    res = np.copy(padded_signal)
    
    for i in range(L, len(padded_signal) - L):
        res[i] = np.mean(padded_signal[i - L:i + L + 1])
    
    # Remove the padded portions
    return res[L:-L]

# Function to read and plot data from multiple files
def plot_data():
    # Open a file dialog to select multiple .txt files
    file_paths = filedialog.askopenfilenames(
        title="Select files",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    
    if not file_paths:
        messagebox.showerror("Error", "No files selected")
        return

    try:
        # Create three separate figures
        plt.figure("Original Spectra", figsize=(10, 8))
        plt.figure("Cropped Spectra", figsize=(10, 8))
        plt.figure("Fit Intersections", figsize=(10, 8))

        for idx, file_path in enumerate(file_paths, start=1):  # Start enumeration from 1
            x_data = []
            y_data = []

            with open(file_path, 'r') as file:
                data_started = False

                for line in file:
                    # Skip everything until "[TRACE DATA]"
                    if "[TRACE DATA]" in line:
                        data_started = True
                        continue

                    # Read the data after "[TRACE DATA]"
                    if data_started:
                        try:
                            # Split the line by comma and strip whitespace
                            x, y = map(float, line.strip().split(','))
                            x_data.append(x)
                            y_data.append(y)
                        except ValueError:
                            # Skip lines that can't be parsed
                            continue
            x_data = np.array(x_data)
            y_data = np.array(y_data)

            L = 10
            y_data = lissage(y_data, L)

            x_lim1 = 1527
            x_lim2 = 1560
            threshold = -75

            mask_x = (x_data >= x_lim1) & (x_data <= x_lim2)
            mask_y = y_data >= threshold

            mask = mask_x & mask_y

            x_cropped = x_data[mask]
            y_cropped = y_data[mask]

            # Plot original data in the first figure
            plt.figure("Original Spectra")
            plt.plot(x_data, y_data, label=f"Spectrum {idx}")

            # Plot cropped data in a second figure
            plt.figure("Cropped Spectra")
            plt.plot(x_cropped, y_cropped, label=f"Spectrum {idx}")

            # Fit lines and find intersection
            n = 200
            max_value = np.max(y_cropped)
            threshold = max_value - 3
            peak_index = np.argmax(y_cropped)

            # Find indices where the spectrum goes below the threshold on the left and right
            left_indices = np.where(y_cropped[:peak_index] <= threshold)[0]
            right_indices = np.where(y_cropped[peak_index:] <= threshold)[0]

            # Define P1 and P2 based on valid indices
            P1 = left_indices[-1] + 1  # Last index below threshold on the left side
            P2 = right_indices[0] + peak_index  # First index below threshold on the right side

            # Get points for regression
            left_x = x_cropped[P1 - n:P1 + n + 1]
            left_y = y_cropped[P1 - n:P1 + n + 1]
            right_x = x_cropped[P2 - n:P2 + n + 1]
            right_y = y_cropped[P2 - n:P2 + n + 1]

            # Perform linear regression for the left side
            model_left = LinearRegression()
            model_left.fit(left_x.reshape(-1, 1), left_y)
            m_left = model_left.coef_[0]
            b_left = model_left.intercept_

            # Perform linear regression for the right side
            model_right = LinearRegression()
            model_right.fit(right_x.reshape(-1, 1), right_y)
            m_right = model_right.coef_[0]
            b_right = model_right.intercept_

            # Create left and right fit functions
            def left_fit(x):
                return m_left * x + b_left

            def right_fit(x):
                return m_right * x + b_right

            # Find intersection of left and right fit lines
            intersection_x = (b_right - b_left) / (m_left - m_right)
            intersection_y = left_fit(intersection_x)

            fit_range = 10

            # Plot the left and right fits in the "Fit Intersections" figure
            left_slope_display_x = np.linspace(intersection_x - fit_range, intersection_x, 100)
            plt.figure("Fit Intersections")
            plt.plot(x_cropped, y_cropped, label=f"Spectrum {idx} - λ = {intersection_x:.2f} nm")
            plt.plot(left_slope_display_x, left_fit(left_slope_display_x), '--', color='blue')
            right_slope_display_x = np.linspace(intersection_x, intersection_x + fit_range, 100)
            plt.plot(right_slope_display_x, right_fit(right_slope_display_x), '--', color='green')
            plt.scatter(intersection_x, intersection_y, color='red')  # Mark intersection with red dot

            # Adding labels, legends, and grid
            plt.xlabel("Wavelength (nm)", fontsize=20)
            plt.ylabel("Power (dBm)", fontsize=20)
            plt.title(f"Fit Intersections", fontsize=25)
            plt.legend(fontsize=20)
            plt.grid(True)
            plt.tick_params(axis='both', which='major', labelsize=20)

        # Show the figures
        plt.show()

    except Exception as e:
        messagebox.showerror("Error", f"Failed to read files: {e}")

# Setup the tkinter GUI
root = tk.Tk()
root.title("Optical Spectrum Analyzer Plotter")

# Create a button to trigger file selection and plotting
button = tk.Button(root, text="Select Files and Plot", command=plot_data)
button.pack(pady=20)

# Set window dimensions
root.geometry("300x100")

# Run the tkinter main loop
root.mainloop()
