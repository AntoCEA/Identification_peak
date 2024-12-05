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

            y_data_final = savgol_filter(y_normalized, window_length=151, polyorder=3, delta=2)


            # Plot original data in the first figure
            plt.figure("Original Spectra")  # Ensure plotting in the original spectrum figure
            plt.plot(x_data, y_data, label=f"Spectrum {idx}")

            # Plot cropped data in a second figure
            plt.figure("Cropped Spectra", figsize=(10, 8))  # Switch to cropped spectrum figure
            plt.plot(x_cropped, y_cropped, label=f"Spectrum {idx}")
            
            plt.figure("Normalized Cropped Spectra", figsize=(10, 8))  # Switch to cropped spectrum figure
            plt.plot(x_cropped, y_normalized, label=f"Spectrum {idx}")

            plt.figure("Normalized Cropped Filtered Spectra", figsize=(10, 8))  # Switch to cropped spectrum figure
            plt.plot(x_cropped, y_data_final, label=f"Spectrum {idx}")

            # Perform Gaussian Fit
            try:
                popt_gauss, _ = curve_fit(
                    gaussian_with_offset,
                    x_cropped,
                    y_data_final,
                    p0=[1, 1545, 5, 0]  # Initial guess
                )
                A0, lambda_max, sigma, offset = popt_gauss

                # Plot Gaussian fit on the cropped data
                plt.figure("Gaussian Fit")
                spectrum_line, = plt.plot(x_cropped, y_data_final)
                plt.plot(
                    x_cropped,
                    gaussian_with_offset(x_cropped, *popt_gauss),
                    '--',
                    color=spectrum_line.get_color(),  # Match the color of the Gaussian fit to the spectrum
                    label=f"Gaussian Fit {idx}: λmax={lambda_max:.2f} nm"
                )

                # Scatter for maximum
                plt.scatter([lambda_max], [A0 + offset], color=spectrum_line.get_color(), marker='o')

                # Correct Longwave Edge (1/e point)
                longwave_edge = lambda_max + sigma  # Longwave edge is one sigma away from lambda_max
                edge_value = A0 / np.e + offset     # 1/e height value
                #plt.scatter([longwave_edge], [edge_value], color=spectrum_line.get_color(), marker='x', label=f'Edge: λ={longwave_edge:.2f} nm')

                plt.legend(fontsize=12)
            except RuntimeError as e:
                print(f"Gaussian fit failed for Spectrum {idx}: {e}")


        # Configure the original spectra plot
        plt.figure("Original Spectra")
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power (dBm)", fontsize=20)
        plt.title("Original Sapphire Bragg Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)

        # Configure the cropped spectra plot
        plt.figure("Cropped Spectra")
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power (dBm)", fontsize=20)
        plt.title("Cropped Sapphire Bragg Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=20)
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)

        # Configure the cropped spectra plot
        plt.figure("Normalized Cropped Spectra")  # Switch to cropped spectrum figure
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power (dBm)", fontsize=20)
        plt.title("Normalized Cropped Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=20)
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)

        # Configure the cropped spectra plot
        plt.figure("Normalized Cropped Filtered Spectra")  # Switch to cropped spectrum figure
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power (dBm)", fontsize=20)
        plt.title("Normalized Cropped Filtered Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=20)
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)

        # Configure the Gaussian fit plot
        plt.figure("Gaussian Fit")
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Normalized Power", fontsize=20)
        plt.title("Gaussian Fit of Cropped Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)

        # Display all figures
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