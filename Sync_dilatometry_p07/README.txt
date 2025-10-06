============================================================
  SYNCHRONIZING DILATOMETRY DATA WITH DIFFRACTION IMAGES
  (P07 BEAMLINE - HEREON)
============================================================

To synchronize your dilatometry data with diffraction images
from experiments performed at the P07 diffraction beamline,
coordinated by the Hereon staff, you can use the executable:

    p07_sync_dila.exe

------------------------------------------------------------
INSTRUCTIONS
------------------------------------------------------------

1. Create a folder containing all numbered metadata files
   from your experiment (including those split due to
   "fast mode", if applicable).

2. Add to the same folder all corresponding .log files
   from your experiment (more than one file may exist if
   "fast mode" acquisition was used).

3. Include the .D5DT file corresponding to your experiment.

4. Run the synchronization tool by either:
      - Double-clicking p07_sync_dila.exe
        OR
      - Executing the script p07_sync_dila.py
        using a suitable Python interpreter.

------------------------------------------------------------
OUTPUT
------------------------------------------------------------

If all files are correctly placed, the program will generate:

 - An image with three graphs showing:
       * Time vs Temperature
       * Time vs Change in Length
       * Time vs Force

 - A synchronized data file named:
       ...sync_file.txt

The sync_file.txt contains the following columns:

   (1) Image name
   (2) Time
   (3) Temperature
   (4) Change in length
   (5) Force

------------------------------------------------------------
NOTE
------------------------------------------------------------

In the ...sync_file.txt, some images may not include time, temperature,
change in length, or force data if they were acquired
before the dilatometer was initialized.
