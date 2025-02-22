import tkinter as tk
from tkinter import PhotoImage
from tkinter import filedialog, messagebox, Toplevel, Text, Button
import subprocess
import threading
import time

#example usage: python GUI_SHINE.py

class GUI_SHINE:
    
    def __init__(self, root):
        """Initialize the GUI."""
        self.root = root
        self.root.title("SHINE (v1.1 version)")
        
        # Set minimum size for the window
        self.root.geometry("1200x1000")  # Width x Height
        self.root.minsize(600, 600)

        # Configure grid to make it expandable
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        
        # Increase default font size for better readability
        self.default_font = ("Arial", 12)
        self.root.option_add("*Font", self.default_font)
        
        
        # Icon as png
        #logo = PhotoImage(file="logo.png")
        #self.root.iconphoto(False, logo)

        

        # === Input Arguments Section (grpinp) ===
        self.frame_input = tk.LabelFrame(root, text="Input Arguments", padx=10, pady=10)
        self.frame_input.grid(row=0, column=0, padx=10, pady=5, sticky="ew")
        
        # Cube and Variance Cube (Browse)
        self.entry_cube    = self.create_input_row_with_browse(self.frame_input, "Cube/Image (Required):", 1)
        self.entry_varcube = self.create_input_row_with_browse(self.frame_input, "Variance Cube/Image (Required):", 2)
        
        # 2D Mask and 3D Mask (Browse)
        self.entry_mask2d     = self.create_input_row_with_browse(self.frame_input, "2D Mask (Optional):", 3)
        self.entry_mask2dpost = self.create_input_row_with_browse(self.frame_input, "2D Mask Post Smoothing (Optional):", 4)
        
        # Small Entries (Cube Extension and Variance Extension)
        self.entry_extcub = self.create_small_input(self.frame_input, "Cube/Image Ext (Default=0):", 0, "0")  
        self.entry_extvar = self.create_small_input(self.frame_input, "Var Cube/Image Ext (Default=0):", 1, "0")  
        
        # Cut the cube and variance if needed
        self.entry_cutzmin = self.create_small_input(self.frame_input, "zmin (Default=None):",2, "None")
        self.entry_cutzmax = self.create_small_input(self.frame_input, "zmax (Default=None):",3, "None")
        self.entry_cutlmin = self.create_small_input(self.frame_input, "lmin (Default=None):",4, "None")
        self.entry_cutlmax = self.create_small_input(self.frame_input, "lmax (Default=None):",5, "None")
        
        # === Extraction Arguments Section (grpext) ===
        self.frame_extraction = tk.LabelFrame(root, text="Extraction Arguments", padx=10, pady=10)
        self.frame_extraction.grid(row=1, column=0, padx=10, pady=5, sticky="ew")

        self.entry_snthresh     = self.create_input_row(self.frame_extraction, "S/N Threshold (Default=2):", 0, "2")
        self.entry_spatsmooth   = self.create_input_row(self.frame_extraction, "Spatial Smoothing (Default=0):", 1, "0")
        self.entry_spatsmoothX  = self.create_input_row(self.frame_extraction, "Smooth X (Optional):", 2)
        self.entry_spatsmoothY  = self.create_input_row(self.frame_extraction, "Smooth Y (Optional):", 3)
        self.entry_specsmooth   = self.create_input_row(self.frame_extraction, "Spectral Smoothing (Default=0):", 4, "0")
        self.entry_connectivity = self.create_input_row(self.frame_extraction, "Connectivity (Default=26):", 5, "26")
        self.entry_maskspatedge = self.create_input_row(self.frame_extraction, "Mask Spatial Edges (Default=20):", 6, "20")
        self.usefftconv         = tk.BooleanVar()
        self.dovarsmooth        = tk.BooleanVar(value=1)
        tk.Checkbutton(self.frame_extraction, text="Use FFT", variable=self.usefftconv).grid(row=7, column=1, sticky="w")
        tk.Checkbutton(self.frame_extraction, text="Smooth Vardata", variable=self.dovarsmooth).grid(row=7, column=2, sticky="w")

        # === Cleaning Arguments Section (grpcln) ===
        self.frame_cleaning = tk.LabelFrame(root, text="Cleaning Arguments", padx=10, pady=10)
        self.frame_cleaning.grid(row=2, column=0, padx=10, pady=5, sticky="ew")

        self.entry_minvox  = self.create_input_row(self.frame_cleaning, "Min Voxels (Default=1):", 0, "1")
        self.entry_mindz   = self.create_input_row(self.frame_cleaning, "Min DZ (Default=3):", 1, "3")
        self.entry_maxdz   = self.create_input_row(self.frame_cleaning, "Max DZ (Default=200):", 2, "200")
        self.entry_minarea = self.create_input_row(self.frame_cleaning, "Min Area (Default=1):", 3, "1")

        # === Output Control Arguments Section (grpout) ===
        self.frame_output = tk.LabelFrame(root, text="Output Control Arguments", padx=10, pady=10)
        self.frame_output.grid(row=4, column=0, padx=10, pady=5, sticky="ew")

        self.var_outdir         = self.create_input_row_with_browse_dir(self.frame_output, "Output directory (Deafault=./):", 0)
        self.var_writelabels    = tk.BooleanVar(value=1)
        self.var_writesmcube    = tk.BooleanVar()
        self.var_writesmvar     = tk.BooleanVar()
        self.var_writesmsnrcube = tk.BooleanVar()
        self.var_writesubcube   = tk.BooleanVar()  
        self.var_writevardata   = tk.BooleanVar()

        tk.Checkbutton(self.frame_output, text="Write Labels", variable=self.var_writelabels).grid(row=1, column=0, sticky="w")
        tk.Checkbutton(self.frame_output, text="Write Smoothed Cube/Image", variable=self.var_writesmcube).grid(row=2, column=0, sticky="w")
        tk.Checkbutton(self.frame_output, text="Write Smoothed Variance", variable=self.var_writesmvar).grid(row=3, column=0, sticky="w")
        tk.Checkbutton(self.frame_output, text="Write Smoothed S/N Cube/Image", variable=self.var_writesmsnrcube).grid(row=4, column=0, sticky="w")
        tk.Checkbutton(self.frame_output, text="Write Subcube", variable=self.var_writesubcube).grid(row=5, column=0, sticky="w")
        tk.Checkbutton(self.frame_output, text="Write Var Data", variable=self.var_writevardata).grid(row=6, column=0, sticky="w")
        
        # === Run Button ===
        self.btn_run = tk.Button(root, text="Run Script", command=self.run_script)
        self.btn_run.grid(row=4, column=0, pady=10)
        
        # Ensure the window is in front and focused when opened
        self.root.after(50, self.bring_to_front)
        

    def bring_to_front(self):
        """Bring the window to the front when opened."""
        self.root.wm_attributes("-topmost", 1)  # Keep window on top
        self.root.after(100, lambda: self.root.wm_attributes("-topmost", 0))  # Remove topmost after 100 ms
        self.root.after(100, self.keep_on_top)
    
    def keep_on_top(self):
        """Keep the window on top until the user interacts with it."""
        self.root.wm_attributes("-topmost", 1)
        self.root.after(100, self.check_focus)
    
    def check_focus(self):
        """Check if the window is focused and stop keeping it on top if focused."""
        if self.root.focus_get() == self.root:
            self.root.wm_attributes("-topmost", 0)  # Stop keeping it on top when focused
        else:
            self.root.after(100, self.keep_on_top)

        


    def create_input_row(self, parent, label, row, default_value=""):
        """Helper function to create a labeled input row."""
        tk.Label(parent, text=label).grid(row=row, column=0, sticky="e", padx=5, pady=2)
        entry = tk.Entry(parent, width=30)
        entry.grid(row=row, column=1, padx=5, pady=2)
        if default_value:
            entry.insert(0, default_value)
        return entry
    
    def create_small_input(self, parent, label, row, default_value=""):
        """Helper function to create a small input field next to an existing row."""
        tk.Label(parent, text=label).grid(row=row, column=3, sticky="e", padx=5, pady=2)  
        entry = tk.Entry(parent, width=5)  
        entry.grid(row=row, column=4, padx=5, pady=2)  
        if default_value:
            entry.insert(0, default_value)
        return entry
        
    
    def create_input_row_with_browse(self, parent, label, row):
        """Helper function to create a labeled input row with a browse button."""
        tk.Label(parent, text=label).grid(row=row, column=0, sticky="e", padx=5, pady=2)
        entry = tk.Entry(parent, width=40)
        entry.grid(row=row, column=1, padx=5, pady=2)
        btn_browse = tk.Button(parent, text="Browse", command=lambda: self.browse_file(entry))
        btn_browse.grid(row=row, column=2, padx=5, pady=2)
        return entry        
            
    def create_input_row_with_browse_dir(self, parent, label, row):
        """Helper function to create a labeled input row with a browse button."""
        tk.Label(parent, text=label).grid(row=row, column=0, sticky="e", padx=5, pady=2)
        entry = tk.Entry(parent, width=40)
        entry.grid(row=row, column=1, padx=5, pady=2)
        btn_browse = tk.Button(parent, text="Browse", command=lambda: self.browse_directory(entry))
        btn_browse.grid(row=row, column=2, padx=5, pady=2)
        return entry
    
    
    def browse_file(self, entry):
        """Open a file dialog and insert the selected file path into the given entry."""
        filepath = filedialog.askopenfilename(title="Select a file")
        if filepath:
            entry.delete(0, tk.END)
            entry.insert(0, filepath)
        
        # Force the GUI to update and regain focus after closing the file dialog
        self.root.after(10, self.root.update_idletasks)  # Ensure window updates
        self.root.after(10, self.root.focus_force)  # Force focus on the window
        
    
    
    def browse_directory(self, entry):
        """Open a directory dialog and insert the selected directory path into the given entry."""
        directory = filedialog.askdirectory(title="Select a directory")
        if directory:
            entry.delete(0, tk.END)
            entry.insert(0, directory)
        
        # Ensure the window updates and gets focus after closing the directory dialog
        self.root.after(10, self.root.update_idletasks)  # Ensure window updates
        self.root.after(10, self.root.focus_force)  # Force focus on the window
        
        
    def show_custom_messagebox(self, title, message):
        """Create a custom message box with adjustable styles."""
        dialog = Toplevel(self.root)
        dialog.title(title)
        dialog.geometry("500x300")  # Set size
        dialog.resizable(False, False)

        # Text widget for the message
        text = Text(dialog, wrap="word", bg="#f9f9f9", fg="#333", font=("Arial", 12), padx=10, pady=10)
        text.insert("1.0", message)
        text.configure(state="disabled")  # Make read-only
        text.pack(expand=True, fill="both", padx=10, pady=10)

        # OK button to close
        ok_button = Button(dialog, text="OK", command=dialog.destroy, bg="#007BFF", fg="white", font=("Arial", 12, "bold"))
        ok_button.pack(pady=10)

        # Modal behavior
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.wait_window()
        
                
        
    def run_script(self):
        

        """Run the external script with the provided parameters."""
        # Collect input arguments
        data       = self.entry_cube.get()
        vardata    = self.entry_varcube.get()
        mask2d     = self.entry_mask2d.get()
        mask2dpost = self.entry_mask2dpost.get()
        extcub     = self.entry_extcub.get()
        extvar     = self.entry_extvar.get()
        zmin       = self.entry_cutzmin.get()
        zmax       = self.entry_cutzmax.get()
        lmin       = self.entry_cutlmin.get()
        lmax       = self.entry_cutlmax.get()

        # Collect extraction arguments
        snthreshold  = self.entry_snthresh.get()
        spatsmooth   = self.entry_spatsmooth.get()
        spatsmoothX  = self.entry_spatsmoothX.get()
        spatsmoothY  = self.entry_spatsmoothY.get()
        specsmooth   = self.entry_specsmooth.get()
        connectivity = self.entry_connectivity.get()
        maskspatedge = self.entry_maskspatedge.get()
        usefftconv   = self.usefftconv.get()
        dovarsmooth  = self.dovarsmooth.get()

        # Collect cleaning arguments
        minvox  = self.entry_minvox.get()
        mindz   = self.entry_mindz.get()
        maxdz   = self.entry_maxdz.get()
        minarea = self.entry_minarea.get()

        # Collect output control arguments
        outdir         = self.var_outdir.get()
        writelabels    = self.var_writelabels.get()
        writesmcube    = self.var_writesmcube.get()
        writesmvar     = self.var_writesmvar.get()
        writesmsnrcube = self.var_writesmsnrcube.get()
        writesubcube   = self.var_writesubcube.get()
        writevardata   = self.var_writevardata.get()
        
        # Construct the command
        command = [
            "python", "SHINE.py",
            data, vardata,
            f"--extdata={extcub}", f"--extvardata={extvar}",
            f"--snthreshold={snthreshold}", f"--spatsmooth={spatsmooth}",
            f"--connectivity={connectivity}", f"--maskspedge={maskspatedge}",
            f"--minvox={minvox}", f"--mindz={mindz}", f"--maxdz={maxdz}", f"--minarea={minarea}",
        ]
        
        # Add optional arguments
        if mask2d:
            command.append(f"--mask2d={mask2d}")
        if mask2dpost:
            command.append(f"--mask2dpost={mask2dpost}")
        if spatsmoothX:
            command.append(f"--spatsmoothX={spatsmoothX}")
        if spatsmoothY:
            command.append(f"--spatsmoothY={spatsmoothY}")
        if specsmooth:
            command.append(f"--specsmooth={specsmooth}")
        if usefftconv:
            command.append(f"--usefftconv={usefftconv}")
        if dovarsmooth:
            command.append(f"--dovarsmooth={dovarsmooth}")

        if zmin and zmax == "None":
            zmin = None
            zmax = None
        if lmin and lmax == "None":
            lmin = None
            lmax = None
        
        #cut
        if zmin and zmax is not None:
            command.append(f"--zmin={zmin} --zmax={zmax}")
        if lmin and lmax is not None:
            command.append(f"--lmin={lmin} --lmax={lmax}")
            
        #save    
        if writelabels:
            command.append("--writelabels")
        if writesmcube:
            command.append("--writesmdata")
        if writesmvar:
            command.append("--writesmvar")
        if writesmsnrcube:
            command.append("--writesmsnr")
        if writesubcube:
            command.append("--writesubcube")        
        if outdir:
            command.append(f"--outdir={outdir}")
        if writevardata:
            command.append("--writevardata")

        # Execute the command
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            self.show_custom_messagebox("Success", result.stdout)
        except subprocess.CalledProcessError as e:
            self.show_custom_messagebox("Error", e.stderr)
        finally:
            # Close the GUI
            self.root.destroy()

def main():

    root = tk.Tk()
    app = GUI_SHINE(root)
    root.mainloop()
    
if __name__ == "__main__":

    main()   
