# Performance of Integrated Relay-RIS in UAV Systems: A Selection Combining Approach

> MATLAB simulation repository

## Overview

This repository implements simulations for the performance analysis of an integrated relay + reconfigurable intelligent surface (RIS) aided unmanned aerial vehicle (UAV) communication system employing selection combining (SC).
The main focus is on metrics like outage probability (OP) and symbol error rate (SER) for different links: direct link, relay-based link, RIS-based link, and UAV-assisted links under Nakagami-m fading and RWP mobility models.

## Project Scope

* Comparison of **Direct**, **DF/AF Relay**, **RIS only**, **UAV Relay/IRS** systems.
* Selection combining among available paths to enhance reliability.
* Analysis under Nakagami-m fading channels, different number of RIS elements, transmit power, noise power, mobility (RWP).
* MATLAB code-based simulation of outage probability and SER.
* Analytical / semi-analytical modelling via Gauss-Laguerre integration weights for performance metrics.

## Repository Structure

```
/                             ← root  
  Direct_OP.m                 ← Outage probability (OP) for direct link  
  Direct_SER.m                ← Symbol error rate (SER) for direct link  
  IRS_OP.m / IRS_OP_new.m     ← OP analysis for RIS-assisted link  
  IRS_SER.m / IRS_SER_new.m   ← SER analysis for RIS-assisted link  
  Relay_OP.m / Relay_OP_new.m ← OP for relay-assisted link  
  Relay_SER.m / Relay_SER_new.m← SER for relay assisted link  
  UAV_OP.m / UAV_OP_ana.m     ← OP for UAV assisted link (simulation / analytical)  
  UAV_SER.m                   ← SER for UAV assisted link  
  OutageNakagamiBPSK.m        ← OP for Nakagami-m fading with BPSK  
  gauss_laguerre_weights.m    ← Function to compute weights for integration  
  RWP.m / RWP_OPTION.m        ← Random Waypoint (mobility) modelling scripts  
  SER.m / SER_IRS.m / SER_Nakagami_m.m / SER_RELAY.m / SER_UAV.m ← General/variant SER scripts  
  TAS.m / TLS.m               ← (If applicable) Transmit antenna selection / transmit-link selection scripts  
```

## Getting Started

### Prerequisites

* MATLAB (version compatible with all used features — e.g., 2018b or later)
* Basic familiarity with wireless communications, fading channels, RIS, UAV mobility.

### Setup & Running Simulations

1. Clone or download this repository:

   ```bash
   git clone https://github.com/Rishikesh24062003/Performance-of-Integrated-Relay-RIS-in-UAV-Systems-A-Selection-Combining-Approach.git
   cd Performance-of-Integrated-Relay-RIS-in-UAV-Systems-A-Selection-Combining-Approach
   ```
2. In MATLAB, add the repository folder to path or set as current folder.
3. Identify the simulation script you want to run. For example:

   * `Direct_OP.m` → runs outage probability for the direct link.
   * `IRS_SER_new.m` → runs SER for RIS-assisted link with updated/new parameters.
   * `UAV_OP_ana.m` → runs analytical outage probability for UAV assisted link.
4. Open the script, modify parameters as needed (transmit power, noise, number of RIS elements, Nakagami-m parameter m, SNR threshold, mobility parameters for RWP).
5. Run the script. Results (curves, tables) will be generated (typically using MATLAB plots). Save figures if required for your report.

### Simulation Parameters

Key parameters you may see in the scripts:

* **Transmit Power** (Pt)
* **Noise Power** (N0)
* **Nakagami-m fading parameters** (m)
* **Number of RIS Elements** (N)
* **SNR threshold** for outage calculation
* **Mobility model parameters** (for RWP: speed range, pause time, simulation area)
* **Selection Combining**: number of available links paths (direct, relay, RIS, UAV) and how SC selects the best.

## Results & Interpretation

* Outage Probability (OP): probability that SNR falls below the threshold → higher number of RIS elements or better relay link should reduce OP.
* Symbol Error Rate (SER): error probability for a given modulation (e.g., BPSK, QPSK) under the considered channel conditions.
* Mobility effects: Random Waypoint (RWP) model affects fading statistics, link quality, and hence OP/SER.
* Selection Combining: shows performance improvement by choosing the best path among multiple available.
* Comparison across system architectures: direct link vs relay vs RIS vs UAV vs combined setups.

## Contribution / Modification

If you plan to extend or adapt this code:

* Add unsupported modulation schemes (e.g., QPSK, 16-QAM) by modifying SER scripts.
* Incorporate more realistic channel models (path loss, shadowing, 3D UAV mobility).
* Extend to multi-UAV or multiple RIS elements/clusters.
* Use different combining techniques (e.g., maximal-ratio combining, equal gain combining) instead of selection combining.
* Automate parameter sweeps and generate summary dashboards.

## Licence & Citation

* The code here is provided “as is”—please use it for academic/research purposes and cite accordingly.
* If you publish results based on this work, please reference the repository and any associated publication.
