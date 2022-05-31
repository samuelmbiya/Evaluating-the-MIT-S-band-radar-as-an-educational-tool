# Evaluating the MIT S-band radar as an educational tool
This code repository contains the code that was written to perform ISAR imaging usng the MIT S-band radar in FMCW mode. This was written for the ISAR demonstrator phase of the final year project.

## Repository information:

- mit_radar_isar_processing.m contains the code used for signal processing and ISAR image processing.

- mit_radar_generate_isar_image.m calls mit_radar_isar_processing files and allows a user to import the recorded radar dataset.

- The datasets folder contains sample test data from the dataset that was used.

- The final_code folder contains a modified version of the code that should be run. See: the update log below

## Instructions to run code:

To successfully execute the ISAR image generating code:

1. Open the final_code folder
2. Include your local file path that points to the dataset folder on line 26 in the mit_radar_generate_isar_image.m script, i.e 

   ```
   All_ISAR_wav_recordings = {[filepath,'1_corner_reflectors_1m_downrange_slow.wav]';[filepath,'2_corner_reflectors_1m_downrange_slow.wav']};
   ```

2. Run the script mit_radar_generate_isar_image.m to generate ISAR images in the recorded radar sample dataset.


## Update log:

On the 31/05/2022 it was found that there were syntax errors in the final code submitted on 30/05/2022 that would prevent the code fom running successfully. The updates to the code are described below:

Script name: mit_radar_generate_isar_image.m

- Line 31 was initially

```
mit_radar_isar_processing(All_ISAR_wav_recordings(ISAR_recording_no)); 
```

was updated to

```
mit_radar_isar_processing(All_ISAR_wav_recordings{ISAR_recording_no}); 
```



Script name: mit_radar_isar_processing.m

- The variable turn_table_speed was missing and added to Line 20 as

```
turn_table_speed = 3.635;
```



- The variable no_pulses_to_consider was missing and was added to Line 19 as

```
no_pulses_to_consider = 16;
```



- Line 26 was initially:

```
All_ISAR_wav_recordings = {'1_corner_reflectors_1m_downrange_slow_FMCW_test2.wav';'2_corner_reflectors_1m_downrange_slow.wav'};
```

and was updated to

```
All_ISAR_wav_recordings = {'1_corner_reflectors_1m_downrange_slow.wav';'2_corner_reflectors_1m_downrange_slow.wav'};
```

This change was performed to reflect the filenames of the original test radar dataset uploaded.
