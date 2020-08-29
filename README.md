# Midi-Voice-Transposer
WIP.  Doesn't really work right now.


Attempts to convert WAV clips of vocals to MIDI files where the end result is recognizable speech, with the goal of improving on existing converters for isolated voice clips.  Very WIP, still includes hardcoded paths to test files etc.

Current status: Seems to read the data correctly and output accurate spectrogram.  Conversion to MIDI sounds bad on any resolution, and playback is glitchy (likely due to the huge size and number of notes in the resulting file).  Suspect some issue with the amplitudes, because they are extremely low on my test file.  May also need to combine multiple contiguous midi notes into a single note, for performance.  Need to create a basic sine wave soundfont to playback as well - at higher resolutions, the attack of the Piano font drowns out the actual note when they're played so quickly
