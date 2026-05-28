## Code Overview

This repository provides the code to reproduce the analysis and figures in the paper:

> **"Computability shapes curiosity: the double helix of scientific convergence"**

### Scripts

- **`analysis.py`**  
  Queries the PubMed database, computes key metrics, and performs Mann‑Kendall trend tests.  
  Outputs:
  - `benchmark_concentration.csv`
  - `benchmark_data_.csv`

- **`visualization.py`**  
  Generates **Figure 2** in the paper using `benchmark_concentration.csv`.

### Requirements

- Python **≥ 3.9**
- Dependencies:
  - `BioPython`
  - `pandas`
  - `numpy`
  - `scipy`
  - `matplotlib`

### How to Run

1. **Set your email** in `analysis.py` (required by Entrez):
   ```python
   Entrez.email = "your-email@example.com"
   ```

2. Run the analysis:
   ```bash
   python analysis.py
   ```

3. Generate the figure:
   ```bash
   python visualization.py
   ```
