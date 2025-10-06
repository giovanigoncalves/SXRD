# 🔄 SYNCHRONIZING DILATOMETRY DATA WITH DIFFRACTION IMAGES  
### *(P07 Beamline – HEREON)*

This guide explains how to **synchronize dilatometry data** with **diffraction images** from experiments performed at the **P07 diffraction beamline**, coordinated by the **HEREON staff**.

---

## ⚙️ Executable
Use the provided executable:


p07_sync_dila.exe


Alternatively, you can run the Python script:


python p07_sync_dila.py


---

## 🧭 Instructions

1. **Create a folder** containing all **numbered metadata files** from your experiment  
   *(including those split due to "fast mode", if applicable).*

2. **Add to the same folder** all corresponding `.log` files from your experiment  
   *(more than one file may exist if "fast mode" acquisition was used).*

3. **Include** the corresponding `.D5DT` file for your experiment.

4. **Run the synchronization tool** by either:
   - 🖱️ Double-clicking **`p07_sync_dila.exe`**  
   - 💻 Executing **`p07_sync_dila.py`** in a suitable Python environment

---

## 📤 Output

If all files are correctly placed, the program will generate:

- 🖼️ **An image** with three graphs showing:
  - ⏱️ *Time vs Temperature*  
  - 📏 *Time vs Change in Length*  
  - ⚙️ *Time vs Force*
- 📄 **A synchronized data file** named:


...sync_file.txt


---

## 📊 sync_file.txt Contents

| Column | Description |
|:------:|--------------|
| (1) | Image name |
| (2) | Time |
| (3) | Temperature |
| (4) | Change in length |
| (5) | Force |

---

## 📝 Note

Some images in the `...sync_file.txt` may **not include** time, temperature, change in length, or force data — this happens if they were acquired **before the dilatometer initialization**.

---

### 🧠 Tip
For best results, make sure:
- File names are consistent and correctly numbered.
- The `.D5DT` file matches the same experiment set.
- All `.log` and metadata files are in the same directory before running the script.

---

✨ **Developed for the P07 Beamline data synchronization workflow (HEREON)**  
📍 *Compatible with standard dilatometry + diffraction experiments*
