# Evaluating the MIT S-band radar as an educational tool
This code repository contains the code that was written to perform Inverse Synthetic Aperture Radar (ISAR) imaging usng the MIT S-band radar in FMCW mode. This was written for the ISAR demonstrator phase of the final year project.

## Repository information:

- mit_radar_isar_processing.m contains the code used for signal processing and ISAR image processing.

- mit_radar_generate_isar_image.m calls mit_radar_isar_processing files and allows a user to import the recorded radar dataset.

- The datasets folder contains sample test data from the dataset that was used.

- The final_code folder contains a modified version of the code that should be run.

- The final_report folder contains the final research report that was written.

## Instructions to run code:

To successfully execute the ISAR image generating code:

1. Open the final_code folder
2. Include your local file path that points to the dataset folder on line 26 in the mit_radar_generate_isar_image.m script, i.e 

   ```
   All_ISAR_wav_recordings = {[filepath,'1_corner_reflectors_1m_downrange_slow.wav]';[filepath,'2_corner_reflectors_1m_downrange_slow.wav']};
   ```

2. Run the script mit_radar_generate_isar_image.m to generate ISAR images in the recorded radar sample dataset.


