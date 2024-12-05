import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import numpy as np

def select_reference_and_plot():
    # Step 1: Select the reference spectrum
    reference_path = filedialog.askopenfilename(
        title="Select Reference Spectrum",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    
    if not reference_path:
        messagebox.showerror("Error", "No reference spectrum selected")
        return

    # Step 2: Select other spectra
    file_paths = filedialog.askopenfilenames(
        title="Select Spectra to Compare",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    
    if not file_paths:
        messagebox.showerror("Error", "No comparison spectra selected")
        return

    try:
        # Load the reference spectrum
        ref_x = []
        ref_y = []
        with open(reference_path, 'r') as ref_file:
            data_started = False
            for line in ref_file:
                if "[TRACE DATA]" in line:
                    data_started = True
                    continue
                if data_started:
                    try:
                        x, y = map(float, line.strip().split(','))
                        ref_x.append(x)
                        ref_y.append(y)
                    except ValueError:
                        continue
        
        ref_x = np.array(ref_x)
        ref_y = np.array(ref_y)

        # Create a figure for the plots
        plt.figure(figsize=(10, 8))

        for idx, file_path in enumerate(file_paths, start=1):
            x_data = []
            y_data = []

            with open(file_path, 'r') as file:
                data_started = False
                for line in file:
                    if "[TRACE DATA]" in line:
                        data_started = True
                        continue
                    if data_started:
                        try:
                            x, y = map(float, line.strip().split(','))
                            x_data.append(x)
                            y_data.append(y)
                        except ValueError:
                            continue
            
            x_data = np.array(x_data)
            y_data = np.array(y_data)

            # Apply the x_data mask to both x_data and y_data
            x_lim1 = 1527
            mask_x = (x_data >= x_lim1)
            x_data = x_data[mask_x]
            y_data = y_data[mask_x]  # Apply the same mask to y_data

            # Interpolate the reference spectrum to match the x_data points
            ref_y_interp = np.interp(x_data, ref_x, ref_y)

            # Calculate the difference
            y_diff = y_data - ref_y_interp

            # Plot the difference
            plt.plot(x_data, y_diff, label=f"Spectrum {idx} - Reference")

        # Configure the plot
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power Difference (dBm)", fontsize=20)
        plt.title("Difference Spectra", fontsize=25)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)
        plt.show()

    except Exception as e:
        messagebox.showerror("Error", f"Failed to process files: {e}")

# Setup the tkinter GUI
root = tk.Tk()
root.title("Spectrum Difference Plotter")

# Create a button to trigger reference selection and plotting
button = tk.Button(root, text="Select Reference and Spectra", command=select_reference_and_plot)
button.pack(pady=20)

# Set window dimensions
root.geometry("300x100")

# Run the tkinter main loop
root.mainloop()
