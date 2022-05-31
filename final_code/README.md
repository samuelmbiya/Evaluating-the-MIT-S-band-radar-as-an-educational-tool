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
