import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt

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

    plt.figure(figsize=(10, 8))  # Create a single figure for all plots
    
    try:
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

            # Plot the data from the current file
            plt.plot(x_data, y_data, label=f"Spectrum {idx}")  # Use index for the label

        # Configure the plot
        plt.xlabel("Wavelength (nm)", fontsize=20)
        plt.ylabel("Power (dBm)", fontsize=20)
        plt.title("Sapphire Bragg Spectrum with mode compler on", fontsize=25)
        plt.legend(loc='best', fontsize=12)  # Automatically adjust legend position
        plt.grid(True)
        plt.tick_params(axis='both', which='major', labelsize=20)
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
