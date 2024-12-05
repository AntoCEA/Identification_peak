import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# Define the Gaussian function with an offset
def gaussian_with_offset(x, A0, lambda_max, sigma, offset):
    return A0 * np.exp(-((x - lambda_max) ** 2) / (2 * sigma ** 2)) + offset

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
        # Create two separate figures for original and cropped spectra
        plt.figure("Original Spectra", figsize=(10, 8))
        
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

            # Define cropping range
            x_lim1 = 1527
            x_lim2 = 1560
            threshold = -75

            mask_x = (x_data >= x_lim1) & (x_data <= x_lim2)
            mask_y = y_data >= threshold

            mask = mask_x & mask_y

            x_cropped = x_data[mask]
            y_cropped = y_data[mask]

            # Normalize y_cropped to the range [0, 1]
            y_min = y_cropped.min()
            y_max = y_cropped.max()
            y_normalized = (y_cropped - y_min) / (y_max - y_min)

            # Define Savitzky-Golay filter parameters
            window_lengths = [121, 131, 151]  # Try different window lengths
            polyorders = [2, 3, 4]  # Try different polynomial orders
            deltas = [1.0, 2.0]  # Try different delta values

            # Create a plot for each combination of parameters
            for window_length in window_lengths:
                for polyorder in polyorders:
                        for delta in deltas:
                            # Apply Savitzky-Golay filter with the current parameters
                            y_filtered = savgol_filter(
                                y_normalized,
                                window_length=window_length,
                                polyorder=polyorder,
                                delta=delta
                            )

                            # Create a subplot to show the filtered data
                            plt.figure(figsize=(10, 8))
                            plt.plot(x_cropped, y_normalized, label="Original", color='blue')
                            plt.plot(x_cropped, y_filtered, label=f"Filtered (window={window_length}, poly={polyorder}, delta={delta})", color='red')
                            plt.xlabel("Wavelength")
                            plt.ylabel("Normalized Intensity")
                            plt.legend()
                            plt.title(f"Spectral Data with Savitzky-Golay Filter\nwindow={window_length}, poly={polyorder}, delta={delta}")
                            plt.grid(True)
                            plt.show()

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")

# Create the tkinter root window
root = tk.Tk()
root.withdraw()  # Hide the main window

# Call the function to plot the data
plot_data()
