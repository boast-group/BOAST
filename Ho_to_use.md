### How to Use the FmpOptBS Toolbox:

1. **Install SPM12:**
   - Ensure that SPM12 is installed. You can download it from the (https://www.fil.ion.ucl.ac.uk/spm/software/download/).

2. **Copy the Toolbox:**
   - Copy the `FmpOptBS` toolbox folder into the `~/SPM/toolbox` directory.

3. **Add the Toolbox to the MATLAB Path:**
   - Add the toolbox directory and all its subfolders to the MATLAB path using the `addpath` function or MATLAB's path manager.

4. **Run SPM for fMRI Analysis:**
   - Launch SPM12 to begin your fMRI analysis.

5. **Start BS Optimization:**
   - Open the SPM batch interface. Under the "SPM" tab, go to "Tools" and select the "BS Optimization" option.

6. **Prepare Your Data:**
   - Add your field map or field gradients, brain mask, and regions of interest (ROIs) using the "Specify" button within the interface.

7. **Execute the Process:**
   - Press the "Play" button in the batch interface to start the optimization process.

8. **Optional: Analyze BS Matrix and ROIs:**
   - If you need to analyze the `BS.matrix` and your selected ROIs, add a pause in the `epi_opt_param_TB` script at line 209.

---

### Modifying Parameters:

You can adjust the fixed protocol parameters in the following scripts:

- **SetDefaultEPIparam.m:**
   - Contains fixed parameters for the EPI protocol. Modify these as needed for your experiment.

- **SetDefaultScannerParam.m:**
   - Defines the default scanner settings. Adjust these to match your scanner's configuration.

- **SetDefaultSimulationParam.m:**
   - Contains default parameters for the simulation environment. You can customize these to fit your specific simulation requirements.

---
