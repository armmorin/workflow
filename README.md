# High-throughput Screening for Low-Temperature Oxide-Ion Conductors

*Predicting oxygen-ion conductivities at different degrees of strain for a crystalline solid through autonomous workflows and Density Functional Theory (DFT).*

---

## 🧭 About This Project

### Real-World Problem

The challenge is to discover and optimize materials that demonstrate **high oxygen-ion conductivity at moderate to lower temperatures** (approximately 400–600 °C) to boost energy efficiency and device longevity. **Ba7Nb4MoO20** and similar hexagonal perovskite-derived materials are promising options owing to their **layered structures**. These structures include palmierite-like and perovskite blocks with inherent oxygen vacancies, which facilitate ion transport.

### Research Motivation

The research aims to **understand and improve oxygen-ion migration mechanisms in these complex materials**. It focuses on how **in-plane strain**, relevant to epitaxial thin-film growth, **impacts ion migration barriers and ionic conductivity**, an aspect that can be strategically managed by choosing different substrates during thin-film fabrication.

Hexagonal perovskites (HPs), such as **Ba7Nb4MoO20, exhibit anisotropic ionic conductivity**. The effects of strain on their crystal structure and oxygen-ion mobility remain underexplored. This study systematically investigates how uniaxial and biaxial strains, both compressive and tensile, influence oxygen migration energy barriers, vacancy formation energies, and local structural distortions, aiming to demonstrate strain engineering as a method for optimizing material performance.

### Significance of the Work

This study shows that in Ba7Nb4MoO20, unlike many other oxygen-ion conductors, **applying compressive in-plane strain *reduces* the migration energy barrier for oxygen ions by as much as 0.14 eV** (~15% decrease), leading to improved ionic conductivity. It identifies specific structural distortions and changes in the coordination environment at the transition states of migrating ions that depend on the strain direction, providing atomic-level insights into the ion migration pathways.

Additionally, the research offers practical advice for the epitaxial growth of HP thin films, indicating that **choosing substrates with certain lattice mismatches (especially inducing compressive strain) can enhance oxygen-ion conductivity**. This strain engineering approach could help design more efficient solid electrolytes for intermediate-temperature fuel cells and sensors.

By emphasizing nanoscale structural control through strain, **the study shifts the focus from solely chemical doping to include mechanical modulation of ion transport**, opening new possibilities for energy devices with better efficiency, stability, and scalability. It combines advanced Density Functional Theory calculations, Nudged Elastic Band techniques, and Continuous Symmetry Measure analyses to develop a comprehensive understanding of the mechano-chemical coupling that governs ionic conductivity in these materials.

In summary, this work addresses the critical energy materials challenge of low-temperature oxide ion conduction by demonstrating strain engineering as an effective method to improve oxygen transport in layered hexagonal perovskite oxides, with important implications for the advancement of next-generation solid-state energy conversion and storage technologies.

## 🔎 Key Insights

- **Compressive in-plane strain consistently lowers the oxygen-ion migration energy barrier** by up to 0.14 eV, resulting in nearly a **15% improvement in ionic conductivity**.

    - This behavior is distinct from many other oxide ion conductors where *tensile strain often enhances conductivity*, highlighting material-specific strain effects.

- The study reveals that **strain induces particular structural distortions and rearrangements in the local coordination environment around migrating oxygen ions**, which facilitate easier ion migration.

    - Continuous Symmetry Measure (CSM) analysis identifies specific polyhedral distortion patterns correlated with reduced migration barriers, providing atomic-scale mechanistic insights.

- The findings suggest that **epitaxial thin-film growth under carefully selected substrate-induced strain can be an effective strategy** to tailor and enhance oxygen-ion conductivity in hexagonal perovskite materials.

- This strain engineering approach offers a **promising avenue beyond chemical doping by mechanically tuning ionic transport properties**, potentially enabling more efficient and stable solid oxide fuel cells and related devices.


## 📝 Methodology

This project employed *state-of-the-art* workflow management tools to investigate oxygen-ion migration in Ba7Nb4MoO20 hexagonal perovskites under strain from first-principles. The methodology can be summarized as follows:

- **Dataset and Systems:**  
  
  Atomic-scale models of Ba7Nb4MoO20 were constructed based on experimentally known crystal structures. The structural information was stored and managed in a SQLite3 database, which was accessed and manipulated using the Atomic Simulation Environment (ASE) toolkit. 
  
  Various strain states were simulated by applying uniaxial (along a or b axes) and biaxial (a and b axes together) in-plane strains, representing both compressive and tensile regimes within ranges relevant to epitaxial thin-film growth (strain levels up to approximately ±3%).

- **Workflow Management**  

  The different stages of the computational workflow, including *structure relaxation*, *vacancy formation energy calculations*, and *Nudged Elastic Band* (NEB) migration barrier calculations, were efficiently orchestrated and managed using PerQueue. This enabled automation, dynamic control, and high-throughput execution of the computational tasks across HPC resources.

- **Calculation Process:**  

  - Oxygen vacancy formation energies were computed to assess the thermodynamic feasibility of vacancy creation under strain. 
  - Oxygen-ion migration barriers along key diffusion pathways were calculated using the Nudged Elastic Band (NEB) method to determine energy barriers for ion hopping between lattice sites.

- **Analytical Tools:**  

  Continuous Symmetry Measure (CSM) analysis quantified distortions in local polyhedral coordination environments around migrating ions, enabling correlation between structural changes and migration barrier variations.


## 📊 Visualizations

![IN PROGRESS!](link_to_key_plot.png)

## 🤔 Interpretation


## 🗂️ Reproducibility

[How could someone else check your results?  
Is your notebook/code available? Is the data public?]: #

**IN PROGRESS**

## 🙏 Acknowledgments / References

-  For more information on this work, you can refer to the preprint version at [Understanding Oxygen-Ion Migration in Hexagonal Perovskites Through Autonomous Workflows and Density Functional Theory](10.26434/chemrxiv-2025-kdhkh)

[-  The database and structures are also publicly available at [Dataset for "Mechanistic Insights into Oxygen-Ion Migration in Hexagonal Perovskites via Autonomous First-Principles Workflows"](10.11583/DTU.29580659)]: #
