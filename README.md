# Autonomous Computational Screening of Hexagonal Perovskites for Enhanced Oxygen-Ion Conductivity

*A high-throughput DFT screening to identify novel oxygen-ion conductors and unravel their atomic-scale migration mechanisms.*

---

## 🧭 About This Project

### Real-World Problem

The development of **sustainable energy technologies** like solid oxide fuel cells (SOFCs) and electrolysis cells (SOECs) **faces a major materials challenge**. These devices need solid electrolytes that can effectively transport oxygen ions, but existing materials such as yttria-stabilized zirconia (YSZ) only show high conductivity at temperatures above 800°C. Running at these high temperatures causes material degradation, compatibility issues with other components, and lowers overall energy efficiency. The main issue is the **absence of materials that provide high oxygen-ion conductivity at lower, more practical temperatures around 400–600°C**.

### Research Motivation

This project was driven by the need to **speed up the discovery of new low-temperature oxygen-ion conducting (OIC) materials**. Recent research has highlighted hexagonal perovskite (HP) structures, such as Ba₇Nb₄MoO₂₀, as promising options due to their inherent oxygen vacancies and adaptable crystal structures. However, a thorough exploration of the extensive chemical space within this structural family was lacking. The main goal was to use autonomous, high-throughput computational workflows to efficiently screen thousands of HP compositions, find new candidates, and **understand the atomic-scale mechanisms that control ion migration**.

### Significance of the Work

This work significantly advances the field of energy materials by systematically screening 5,400 hexagonal perovskite compositions through an entirely autonomous DFT-based workflow. It successfully identified **29 unique and promising candidate materials** for low-temperature oxygen-ion conduction, many of which have not yet been synthesized experimentally.

In addition to discovering new materials, the study offers valuable insights into the physics of ion transport. It clearly categorizes oxygen-ion migration into two separate pathways: a straightforward **ion-hopping** mechanism and a more intricate **cooperative migration** mechanism. Importantly, the research shows that the cooperative mechanism, characterized by the collective rotation of flexible polyhedral units, consistently results in lower energy barriers and improved ion mobility. These insights provide a rational framework for guiding future experimental efforts, emphasizing the design of materials that promote these low-energy cooperative pathways.

## 🔎 Key Insights

- The high-throughput screening effectively reduced a chemical space of 5,400 structures to **29 novel and promising candidate compositions** for low-temperature oxygen-ion conduction.

- Oxygen-ion migration mainly occurs via two mechanisms: a straightforward **ion-hopping** process and a **cooperative mechanism** involving the collective rotation of tetrahedral units.

  - The **cooperative migration** is consistently linked to lower activation barriers because rotating flexible polyhedral units helps minimize high-energy lattice strain during ion transport.

- Analysis using Continuous Symmetry Measures (CSM), the study found a direct link between local coordination, polyhedral distortion, and the migration mechanism, providing a quantitative measure of ionic mobility.

- Key chemical trends that promote the rotational dynamics necessary for the low-energy cooperative mechanism were detected in this study.

## 📝 Methodology

This project was executed using a multi-stage, autonomous computational workflow designed to efficiently screen a large materials database for promising oxygen-ion conductors.

- **Dataset and Systems:**  
  
  The study began with a chemical space of 5,400 unique compositions based on the A₇B₄B'O₂₀ hexagonal perovskite prototype. All structural data, calculations, and results were stored and managed in a **SQLite3 database**, accessed and manipulated using the **Atomic Simulation Environment (ASE)**
  
- **Workflow Management**  

  The entire high-throughput screening process—from structure generation to final analysis—was orchestrated using **PerQueue**, a dynamic Python-based workflow manager. PerQueue enabled the automation of thousands of parallel calculations on high-performance computing (HPC) clusters, managing dependencies, and dynamically filtering candidates based on calculated properties.

- **Calculation Process:**  

  First-principles calculations were performed using **Density Functional Theory (DFT)** as implemented in the VASP code. The workflow consisted of several filtering steps:
  
    1.  **Initial Screening:** Structures were filtered based on chemical and geometric rules (e.g., charge neutrality, ionic radii).
    2.  **Stability Analysis:** Thermodynamic stability was assessed by calculating the energy above the convex hull (E_hull < 100 meV/atom).
    3.  **Electronic Properties**: A bandgap filter (> 1.0 eV) was applied to select for electronically insulating materials.
    4.  **Kinetic Analysis:** For the remaining candidates, the **Nudged Elastic Band (NEB)** method was used to calculate the oxygen-ion migration energy barrier (E_A < 1.1 eV)

-  **Analytical Tools**

  Continuous Symmetry Measure (CSM) analysis was employed to quantify the distortion of the coordination polyhedra around the migrating ion at each step of the NEB path. This allowed for the classification of migration mechanisms and provided a direct link between atomic structure and ionic mobility.

## 📊 Visualizations

**Workflow Funnel**
![A diagram illustrating the multi-step screening process, showing the number of candidate materials that pass or fail at each stage of the computational workflow (e.g., Stability Check, Bandgap, NEB).](https://github.com/armmorin/workflow/blob/main/Sankey_wf.pdf)

**Migration Mechanism Schematics**
![Side-by-side atomic visualizations of the "ion-hopping" and "cooperative" migration pathways. These diagrams would show the evolution of the local polyhedra and the trajectory of the migrating oxygen ion, clearly distinguishing the rotational motion in the cooperative mechanism from the simple hop.](https://github.com/armmorin/workflow/blob/main/ces_ion_hop_coop.pdf)

## 🤔 Interpretation

This screening project's results strongly indicate that tailoring the local geometry and dynamic flexibility of the crystal lattice is more effective for improving ionic conductivity than merely adjusting composition. The main insight is that materials supporting cooperative migration—where polyhedral units rotate to facilitate ion movement—perform better. This offers a clear design principle for experimentalists: prioritize compositions with flexible tetrahedral units and ample free volume to enable these low-energy rotational modes.

Additionally, this research confirms the value of autonomous high-throughput workflows as essential tools in modern materials discovery. By automating complex, repetitive computational screening tasks, we can navigate extensive chemical spaces and identify fundamental structure-property relationships that would be difficult to uncover with traditional experimental or computational approaches.

## 🗂️ Reproducibility

- **Code Availability**: The core of the project is built on open-source tools. The workflow automation was managed by PerQueue, and all structure/data handling was done with the Atomic Simulation Environment (ASE). The specific scripts for generating, filtering, and analyzing the data in this project can be made available in a repository.

- **Data Accessibility**: The full dataset, including the initial structures, DFT calculation results, and the final 29 candidate materials, will be publicly available for verification and further analysis. A dedicated DOI for the dataset can be provided.

- **Environment and Dependencies**: The project was run in a Linux-based HPC environment. Key software dependencies include:

  -   Python (>3.10)
  -   PerQueue (v0.3.0a0)
  -   ASE (v3.22.1)
  -   Pymatgen
  -   VASP (v6.4.2) for DFT calculations

- **Instructions to Reproduce Key Results**: The workflow can be reproduced by installing the required dependencies, configuring the `project_config.cfg` file with the appropriate calculation parameters, and launching the main workflow script via PerQueue. Subsequent analysis scripts can then be run on the output database to regenerate the figures and key findings.

- **Expected Outcomes**: A successful run of the workflow will produce a populated SQLite3 database containing the calculated properties for all structures and a final filtered list of the 29 candidate materials with low migration barriers.

## 🙏 References / Data Availability

-  The primary findings of this work are detailed in the scientific article: [Understanding Oxygen-Ion Migration in Hexagonal Perovskites Through Autonomous Workflows and Density Functional Theory.](https://doi.org/10.26434/chemrxiv-2025-kdhkh)

-  The workflow manager used is described in: [PerQueue: Dynamic Workflow Manager for Materials Discovery.](https://doi.org/10.1039/D4DD00134F)

-  The complete dataset, including structures and calculation outputs, will soon be available.
