import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



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
            y_data1 = np.array(y_data)

            L = 10
            y_data = lissage(y_data1, L)

            x_lim1 = 1527
            x_lim2 = 1560
            threshold = -75

            mask_x = (x_data >= x_lim1) & (x_data <= x_lim2)
            mask_y = y_data >= threshold

            mask = mask_x & mask_y

            x_cropped = x_data[mask]
            y_cropped = y_data[mask]

            # Plot original data in the first figure
            plt.figure("Original Spectra")  # Ensure plotting in the original spectrum figure
            plt.plot(x_data, y_data, label=f"Spectrum {idx}")

            # Plot cropped data in a second figure
            plt.figure("Cropped Spectra", figsize=(10, 8))  # Switch to cropped spectrum figure
            plt.plot(x_cropped, y_cropped, label=f"Spectrum {idx}")

            # Polynomial Fit with Order 4
            initial_guess_poly = [0, 0, 0, 0, y_cropped.min()]  # Initial guess for quartic coefficients
            popt_poly, _ = curve_fit(
                lambda x, a, b, c, d, e: a * x**4 + b * x**3 + c * x**2 + d * x + e,
                x_cropped, y_cropped,
                p0=initial_guess_poly,
                maxfev=100000
            )

            # Extract polynomial coefficients
            a, b, c, d, e = popt_poly

            # Calculate the maximum of the polynomial
            # Solve the derivative = 0 for a quartic polynomial
            roots = np.roots([4 * a, 3 * b, 2 * c, d])  # Quartic derivative coefficients
            real_roots = [r for r in roots if np.isreal(r) and x_cropped.min() <= r <= x_cropped.max()]
            lambda_max_poly_fit = real_roots[0].real if real_roots else np.nan

            # Define the quartic polynomial function for plotting
            def quartic_polynomial(x, a, b, c, d, e):
                return a * x**4 + b * x**3 + c * x**2 + d * x + e
            
            colors = plt.cm.plasma(np.linspace(0, 1, len(file_paths)))  # Use colormap for distinct colors


            # Plot Polynomial Fit
            plt.figure("Polynomial Fit (Order 4)")
            plt.plot(x_cropped, y_cropped, color=colors[idx-1])  # Spectrum color
            plt.plot(x_cropped, quartic_polynomial(x_cropped, *popt_poly), '--', 
                    label=f"Poly Fit {idx}: λmax={lambda_max_poly_fit:.2f} nm", color=colors[idx-1])  # Fit color
            if real_roots:
                plt.plot(lambda_max_poly_fit, quartic_polynomial(lambda_max_poly_fit, *popt_poly), 
                        'o', color=colors[idx-1], markersize=8)  # Maximum point color
            plt.legend(fontsize=12)
            plt.xlabel("Wavelength (nm)", fontsize=20)
            plt.ylabel("Power (dBm)", fontsize=20)
            plt.legend(loc='best', fontsize=20)
            plt.tick_params(axis='both', which='major', labelsize=20)
            plt.title("Polynomial Fit", fontsize=25)
            plt.grid(True)

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

        # Display both figures
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
