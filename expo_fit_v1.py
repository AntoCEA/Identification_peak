import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# Define the piecewise fitting function
def piecewise_exponential_fit(wavelength, I0, A, omega, lambda_B0):
    """
    Piecewise function to model FBG reflection peaks.
    The function behaves differently before and after lambda_B0.
    """
    # Function for wavelengths after lambda_B0
    def after_lambda_B0(wavelength):
        exponent = (wavelength - lambda_B0) / omega
        return I0 + A * np.exp(-np.exp(exponent) - exponent + 1)
    
    # Function for wavelengths before lambda_B0 (mirrored version)
    def before_lambda_B0(wavelength):
        exponent = (lambda_B0 - wavelength) / omega
        return I0 + A * np.exp(-np.exp(exponent) - exponent + 1)
    
    # Apply the appropriate function depending on the wavelength
    return np.piecewise(wavelength, 
                        [wavelength < lambda_B0, wavelength >= lambda_B0], 
                        [before_lambda_B0, after_lambda_B0])

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
        # Create separate figures for original and cropped spectra
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

            # Plot cropped, normalized, and filtered spectra
            plt.figure("Normalized Cropped Filtered Spectra", figsize=(10, 8))
            plt.plot(x_cropped, y_data_final, label=f"Spectrum {idx}")



            # Perform Piecewise Exponential Fit
            try:
                initial_guess = [0.0, 1.0, 15, 1547]  # Initial parameter estimates
                popt, _ = curve_fit(
                    piecewise_exponential_fit, 
                    x_cropped, 
                    y_data_final, 
                    p0=initial_guess, 
                )
                I0, A, omega, lambda_B0 = popt

                # Plot the fitted curve
                plt.figure("Piecewise Exponential Fit")
                spectrum_line, = plt.plot(x_cropped, y_data_final, label=f"Data {idx}")
                plt.plot(
                    x_cropped,
                    piecewise_exponential_fit(x_cropped, *popt),
                    '--',
                    color=spectrum_line.get_color(),
                    label=f"Fit {idx}: λB0={lambda_B0:.2f} nm, ω={omega:.2f}"
                )
                plt.plot(x_cropped, piecewise_exponential_fit(x_cropped, *initial_guess), label="Initial Guess", linestyle='dotted')
                print(f"Fitted parameters: I0={I0}, A={A}, ω={omega}, λB0={lambda_B0}")
                plt.legend(fontsize=12)
                plt.xlabel("Wavelength (nm)", fontsize=14)
                plt.ylabel("Normalized Intensity", fontsize=14)
                plt.title("Piecewise Exponential Fit", fontsize=16)
                plt.grid(True)
            except RuntimeError as e:
                print(f"Piecewise exponential fit failed for Spectrum {idx}: {e}")

        # Configure the plots
        plt.figure("Original Spectra")
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power (dBm)", fontsize=20)
        plt.title("Original Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)

        plt.figure("Normalized Cropped Filtered Spectra")
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Normalized Power", fontsize=20)
        plt.title("Normalized Cropped and Filtered Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)

        plt.figure("Piecewise Exponential Fit")
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Normalized Power", fontsize=20)
        plt.title("Piecewise Exponential Fit of Cropped Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)

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
