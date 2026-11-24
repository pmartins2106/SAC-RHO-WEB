# SAC-RHO 🧬

A Streamlit web application for the automatic analysis of seed amplification assays (SAA) and protein aggregation kinetics.

## Features

- **Automated Curve Analysis:** Estimates baseline ($F_0$), endpoint ($F_\infty$), and initial kinetic parameters across multiple protein aggregation curves.
- **SAC & RHO Determination:** Calculates the global Seeding Activity Coefficient (SAC) and individual $\rho$ (RHO) values with propagated standard errors.
- **Customizable Kinetic Parameters:** Configure primary nucleation rates ($k_b$), SAC normalization, reference endpoints, and minimum datapoint thresholds ($N_{\min}$).
- **Interactive Visualization:** Generates dark-mode plots comparing experimental data against fitted kinetic curves.
- **Built-in Template & Example:** Includes a downloadable `.ods` spreadsheet template and a one-click prefilled example dataset.


## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/pmartins2106/SAC-RHO-WEB.git
   cd SAC-RHO-WEB
   ```

2. **Create and activate a virtual environment (recommended):**
   ```bash
   python -m venv .venv
   # On Windows:
   .venv\Scripts\activate
   # On macOS/Linux:
   source .venv/bin/activate
   ```

3. **Install the required dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Start the Streamlit web application:
   ```bash
   streamlit run home.py
   ```

2. Open the local URL shown in your terminal (typically `http://localhost:8501`).

3. **Run an analysis:**
   - Check **Run a prefilled example** to test the tool immediately, **or**
   - Check **Download template**, fill in your time-series fluorescence data in `.ods` format, and upload it via the file uploader.
   - Adjust optional parameters ($k_b$, SAC normalization, reference endpoint, $N_{\min}$) in the left sidebar.

## Input Data Format

- Input files must be saved in OpenDocument Spreadsheet (`.ods`) format matching `datasets/SAC-RHO_template.ods`.
- Column 1 represents **Time**, and subsequent columns represent individual aggregation progress curves (up to 100 curves and 1,000 time points per run).
