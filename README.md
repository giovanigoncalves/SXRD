# 🧩 SXRD

**SXRD** is a project developed to assist in the **analysis of synchrotron X-ray diffraction (SXRD)** data obtained from experiments performed at the **P07 diffraction beamline**, coordinated by **Hereon**.

---

## 📖 Table of Contents
- [Available Subprojects](#-available-subprojects)
  - [Sync_dilatometry_p07](./Sync_dilatometry_p07/README.md)
  - [dislocation_density](./dislocation_density/README.md)
- [About](#-about)

---

## 📦 Available Subprojects

### 1. 🔗 `Sync_dilatometry_p07`

Synchronizes **dilatometry** and **diffraction** data obtained during dilatometry experiments, whether or not a **deformation module** was used.  
This tool allows users to align experimental measurements in time, making it easier to interpret correlated thermal–mechanical–structural data.

📘 **Documentation:** See the [`Sync_dilatometry_p07` README](./Sync_dilatometry_p07/README.md) for detailed usage instructions.

---

### 2. ⚙️ `dislocation_density`

Calculates the **dislocation density** from given **Full Width at Half Maximum (FWHM)** values.  
To use this module, the user must first **refine the diffraction data** and obtain **FWHM results** (e.g., from peak fitting or Rietveld analysis).

We recommend using the open-source software **[Pydidas](https://github.com/hereon-github/pydidas)** (also developed by Hereon) for this purpose.

📘 **Documentation:** See the [`dislocation_density` README](./dislocation_density/README.md) for further details.

---

## 🧠 About

This repository is part of a broader effort to provide **open tools for synchrotron-based materials research**, enabling efficient data handling, synchronization, and microstructural analysis.

---

### 🧑‍💻 Maintainers
Project coordinated and maintained by **Hereon – P07 Beamline Team**.


