using Melanchall.DryWetMidi.Common;
using Melanchall.DryWetMidi.Core;
using Melanchall.DryWetMidi.Interaction;
using Melanchall.DryWetMidi.MusicTheory;

using NAudio.Dsp;
using NAudio.Wave;

using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Text;

using Note = Melanchall.DryWetMidi.Interaction.Note;

namespace Midi_Transposer_Console
{
    /// <summary>
    /// Handles conversion of sound files to MIDI, each instance of a Converter has its own settings etc
    /// </summary>
    public class Converter
    {
        private readonly float[] MidiFrequencyMap;
        private float[] peaks;

        public Converter()
        {
            MidiFrequencyMap = GetMidiFrequencyMap();
            Console.WriteLine();
        }

        public void ConvertWavToMidi(string inputFile)
        {
            // 2048,4096, or 8192
            uint windowSize = 16384;
            float noteThreshold = 0.01f; // Amplitudes above this get translated to a midi note
            // Note that we can translate this amplitude to a velocity though, but we still don't want 10billion unheard notes

            List<float> audioContentList = new List<float>(); // TOOD: Make this an array if we can calc the size
            uint sampleRate = 0;

            float volMultiplier = 1;
            TimeSpan wavDuration = TimeSpan.Zero;
            using (var reader = new AudioFileReader(inputFile))
            {
                sampleRate = (uint)reader.WaveFormat.SampleRate;
                wavDuration = reader.TotalTime;

                float[] buffer = new float[reader.WaveFormat.SampleRate];

                int numSamplesRead;
                float max = 0;
                do
                {
                    numSamplesRead = reader.Read(buffer, 0, buffer.Length);
                    // Find the peak so we can normalize
                    for (int n = 0; n < numSamplesRead; n++)
                    {
                        var abs = Math.Abs(buffer[n]);
                        if (abs > max) max = abs;
                    }
                    audioContentList.AddRange(buffer); // And add it to buffer
                } while (numSamplesRead > 0);


                if (max == 0)
                    throw new InvalidOperationException("No audio data found");

                volMultiplier = 1.0f / max;
                // Multiply each value in the array by this to normalize
            }
            Console.WriteLine($"Normalizing - Multiply volume values by {volMultiplier}");
            float[] audioContent = audioContentList.Select(a => a * volMultiplier).ToArray();
            var spectrumList = Fft(windowSize, audioContent, sampleRate);


            // Just to confirm the data is there and somewhat the right data, let's draw it
            // On a spectrograph, X axis is time, Y axis is frequency, and color is amplitude
            int imageSize = 1024;

            Bitmap b = new Bitmap(imageSize, imageSize);
            int timesteps = spectrumList.Count();

            float xAxisRatio = (float)imageSize / timesteps;
            float yAxisRatio = (float)(imageSize / MidiFrequencyMap.Length);
            double minDb = 0;
            double maxDb = double.MinValue;
            float maxAmp = 0; // TODO: Pick one and remove the other
            float maxAdb = float.MinValue;

            using (Graphics g = Graphics.FromImage(b))
            {
                g.FillRectangle(Brushes.Black, 0, 0, imageSize, imageSize);
                for (int i = 0; i < timesteps; i++)
                {
                    int count = 0;
                    foreach (float j in MidiFrequencyMap)
                    {
                        var amplitude = spectrumList[i][j];
                        var db = 20 * Math.Log10(amplitude);
                        
                        if (db < minDb)
                            minDb = db;
                        if (db > maxDb)
                            maxDb = db;
                        if (amplitude > maxAmp)
                            maxAmp = amplitude;
                        var adb = AWeight(amplitude);
                        if (adb > maxAdb)
                            maxAdb = (float)adb;
                        //amplitude *= 25; // Scale it up so you can see anything for basic display
                        // While running we're calculating the amp normalization for the MIDI to use
                        float x = i * xAxisRatio;
                        float y = imageSize - (float)(count++ * yAxisRatio);
                        //if (amplitude > noteThreshold)
                        //    Console.WriteLine(amplitude);
                        //int color = (int)(amplitude * 255);
                        int color = (int)Math.Clamp(255 + db * 3, 0, 255);

                        // TODO: Change this all to lockbits stuff if it ever becomes a perf issue.  Or just remove it.
                        if (color > 0)
                            g.DrawLine(new Pen(Color.FromArgb(color, color, color), 1), x, y, x + 1, y + 1);

                        if(i == 0)
                        {
                            // Write the frequency at the x,y
                            g.DrawString(j.ToString(), new Font(FontFamily.GenericSansSerif, 8f), Brushes.White, x, y);
                        }
                    }
                }
                g.Save();
            }
            // Annoyingly, we have to loop it again to get the max db values after applying this amplitude adjustment
            // ... Or, just don't apply it.  
            /*
             * 
             * var amplitude = spectrumList[i][j];
                        var db = 20 * (float)Math.Log10(amplitude);
                        if (AWeight(db) < minDb)
                            minDb = db;
            */



            float dbMult = 1.0f / (float)minDb;
            float ampMult = 1.0f / maxAmp;
            Console.WriteLine($"Db Mult: {dbMult}");

            b.Save("test.png", ImageFormat.Png);

            Console.WriteLine("PNG generated, creating midi");
            var midiFile = new MidiFile();
            TempoMap tempoMap = midiFile.GetTempoMap();

            var trackChunk = new TrackChunk();

            // Alright, current goal: Combine contiguous notes.
            // We should store a Note[128] of all notes from the last frame
            // If we're going to add a note and one already exists, just extend the duration
            // If any notes exist but we aren't adding them, finalize adding them to the notes list and remove them from the arr

            Note[] noteBuffer = new Melanchall.DryWetMidi.Interaction.Note[128];


            using (var notesManager = trackChunk.ManageNotes())
            {
                NotesCollection notes = notesManager.Notes;
                for (int i = 0; i < timesteps; i++)
                {
                    uint count = 0; // We need this because j is no longer constrained 0-128
                    foreach (float j in MidiFrequencyMap)
                    {
                        var amplitude = Math.Clamp(spectrumList[i][j], 0, 1);
                        var decibels = 20 * Math.Log10(amplitude); // Apparently
                        var aDecibels = AWeight(amplitude);
                        // This will be negative


                        // Something's not right here
                        // It seems high notes are a lot quieter than low ones when they shouldn't be
                        // Which sounds like a log scale human hearing thing
                        // ... Or just the piano soundfont being weird

                        // I tried using A-Weighting, but it came out way worse, and that really shouldn't be necessary
                        // Something else must be wrong.

                        // Mostly I think I need to convert amplitude to db and then to velocity

                        //var velocity = (float)Math.Clamp(Math.Log10(amplitude+1) * 0.33 * 127, 0, 127);
                        //var velocity = Math.Pow(2, amplitude * 7);

                        // So I've tried a million different ways, this is the only one that sounds half decent
                        //var velocity = (float)Math.Clamp((Math.Log10(amplitude) + 1) * 127, 0, 127); // * ampMult?  * dbMult?  I think just an arbitrary scaling here is fine
                        // As far as I can tell, amplitude is exponential (2^x) from velocity according to some papers
                        // So that 2^velocity should be amplitude.  So Log2(amplitude) should be velocity
                        // Using Log10 because it simplifies things and should be basically the same thing (sounded better than Log2)


                        // But, maybe, are these amplitude values related to the peaks from the first pass?
                        // Like will all amplitudes for a timestep add up to 1.0
                        // Or will they add up to the value of the peak?  



                        //var velocity = amplitude * ampMult * 127;
                        //var velocity = (1 - ((decibels - maxDb) / minDb)) * 127;

                        // There are weird situations where, basically if too many notes match, the whole thing comes out at 2kb and doesn't send the NoteOff
                        // Basically with a NoteOn for every key for the entire duration
                        // I guess.  I don't know.  It just breaks sometimes depending on what condition I put here.  It's so odd. 


                        //var velocity = Math.Clamp(127 + (20 * (Math.Log2(amplitude)/3) * 1.5), 0, 127);
                        var velocity = Math.Clamp(127 + (decibels*2 - maxDb), 0, 127);

                        //if (amplitude > noteThreshold)
                        //if (aDecibels > decibelThreshold)
                        if ((SevenBitNumber)velocity > 20)
                        //if (decibels > maxDb * 4)
                        {
                            //var velocity = (1 - ((decibels - maxDb) / (maxDb))) * 127;

                            Console.WriteLine(velocity);
                            // Convert i into an absolute time...
                            // Which means, taking the length of the wav/i

                            double durationMs = (wavDuration / timesteps).TotalMilliseconds;
                            double positionMs = i * durationMs;

                            if (noteBuffer[count] != null)
                            {
                                
                                // There was a contiguous note before this
                                var prevNote = noteBuffer[count];
                                // Edit the note, adding durationMs to the length
                                prevNote.Length += LengthConverter.ConvertFrom(
                                    new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)durationMs),
                                    0,
                                    tempoMap);
                                // And... I guess, we should make the Velocity the average of them all
                                // This is probably the part that's going to cause problems

                                prevNote.Velocity = (SevenBitNumber)((prevNote.Velocity + velocity) / 2);
                                

                                /*
                                // Temporarily, save the old note and store a new one
                                notes.Add(noteBuffer[count]);
                                var noteNum = (SevenBitNumber)GetMidiNoteNumber(j);
                                //Console.WriteLine($"Freq {j} becomes {noteNum}");
                                var note = new Note(
                            noteNum,
                            LengthConverter.ConvertFrom(
                                new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)durationMs),
                                0,
                                tempoMap),
                            LengthConverter.ConvertFrom(
                                new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)positionMs),
                                0,
                                tempoMap));
                                note.Velocity = (SevenBitNumber)velocity;
                                noteBuffer[count] = note;
                                */
                            }
                            else
                            {
                                var noteNum = (SevenBitNumber)GetMidiNoteNumber(j);
                                //Console.WriteLine($"Freq {j} becomes {noteNum}");
                                var note = new Note(
                            noteNum,
                            LengthConverter.ConvertFrom(
                                new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)durationMs),
                                0,
                                tempoMap),
                            LengthConverter.ConvertFrom(
                                new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)positionMs),
                                0,
                                tempoMap));
                                note.Velocity = (SevenBitNumber)velocity;
                                noteBuffer[count] = note;
                            }

                        }
                        else
                        {
                            // If the note was off, we need to finalize any previous notes
                            if (noteBuffer[count] != null)
                            {
                                var prevNote = noteBuffer[count];
                                notes.Add(prevNote);
                                noteBuffer[count] = null;
                            }
                        }
                        count++;
                    }
                }
                // We need to finalize anything left in noteBuffer
                foreach (var n in noteBuffer)
                {
                    if (n != null)
                        notes.Add(n);
                }
            }

            midiFile.Chunks.Add(trackChunk);
            midiFile.Write("Vocal_Conversion.mid", true);
            // TODO: Use Polyphone (standalone app) to create a soundfont that is a plain tone
            // And either package this with the midi somehow, or just, use it for playback in midieditor or wherever
        }

        public double AWeight(double f)
        {
            return 20 * Math.Log10(RA(f)) - 20 * Math.Log10(RA(1000));
        }

        double RA(double f)
        {
            double a1 = Math.Pow(f, 2) + Math.Pow(20.6, 2);
            double a2 = Math.Sqrt((Math.Pow(f, 2) + Math.Pow(107.7, 2)) * (Math.Pow(f, 2) + Math.Pow(737.9, 2)));
            double a3 = Math.Pow(f, 2) + Math.Pow(12200, 2);
            return (Math.Pow(12200, 2) * Math.Pow(f, 4)) / (a1 * a2 * a3);
        }

        private int GetMidiNoteNumber(float frequency)
        {
            // Realistically we need to binary search this, this could be slow af otherwise
            // But existing binary search algos won't really work
            // Because it probably won't be an exact match
            // I guess we need to find the one with the lowest difference from our frequency

            // and I can't really figure out how to binary search that atm even though that'd be way better...
            // Well, this should work
            var index = BinarySearch(MidiFrequencyMap, frequency);
            return index;
        }

        public static int BinarySearch(float[] a, float item)
        {
            int first = 0;
            int last = a.Length - 1;
            int mid = 0;
            do
            {
                mid = first + (last - first) / 2;
                if (item > a[mid])
                    first = mid + 1;
                else
                    last = mid - 1;
                if (a[mid] == item)
                    return mid;
            } while (first <= last);
            if (mid == a.Length - 1 || Math.Abs(item - a[mid + 1]) > Math.Abs(item - a[mid]))
                return mid; // It could be either this or mid+1
            else // So we see which one has the lowest difference and return that one.
                return mid + 1;
        }

        private float[] GetMidiFrequencyMap()
        {
            // Calculates a MIDI frequency map such that the array is 128 floats long
            // Each index is the MIDI note number, and each value is the frequency for that note number
            // We just generate this once on instantiation and use that


            float[] result = new float[128];
            for (int i = 0; i < 128; i++)
            {
                result[i] = 8.1758f * (float)(Math.Pow(2, i / 12f)); // Internet said this is how to get them
            }
            return result;
        }

        private FrequencySpectrum ToFrequencySpectrum(Complex[] timeDomain, uint sampleRate)
        {
            int m = (int)Math.Log(timeDomain.Length, 2); // So, either timeDomain.Length needs to be a power of 2
            // Or this needs to be changed to take windowSize instead.  Not sure


            //true = forward fft
            FastFourierTransform.FFT(true, m, timeDomain);
            return new FrequencySpectrum(timeDomain, sampleRate, (uint)timeDomain.Length);
        }

        public IList<FrequencySpectrum> Fft(uint windowSize, float[] audioContent, uint sampleRate)
        {
            IList<Complex[]> timeDomainChunks = this.SplitInChunks(audioContent, windowSize);
            return timeDomainChunks.Select(x => ToFrequencySpectrum(x, sampleRate)).ToList();
        }

        private IList<Complex[]> SplitInChunks(float[] audioContent, uint chunkSize)
        {
            IList<Complex[]> splittedContent = new List<Complex[]>();


            for (uint k = 0; k < audioContent.Length; k += chunkSize)
            {
                long size = k + chunkSize < audioContent.Length ? chunkSize : audioContent.Length - k;
                Complex[] chunk = new Complex[size];

                // This is harder than it seems.
                // I need each chunk to contain the last half of the previous chunk, as well as its own data 
                // And then also half of the following chunk
                // Which means our actual chunks would end up being twice the size that we specify

                // So this can actually be accomplished by just subtracting chunk.Length/4 where possible
                
                // If problematic, can be reverted by just removing that and the if.  Seems okay though.

                for (int i = 0; i < chunk.Length; i++) 
                {
                    // TODO: Figure out the difference in the windows and pick the right one here
                    // Also now let's try adding offsets so it takes half of the previous sample
                    // If we can sub half a chunksize, do it
                    // ... I wonder if that's what the window is doing... 
                    // Anyway, that shouldn't happen here. 
                    var audioIndex = k + i - chunk.Length / 4;
                    if (audioIndex >= 0)
                        chunk[i].X = (float)(audioContent[audioIndex] * FastFourierTransform.BlackmannHarrisWindow((int)audioIndex, (int)chunkSize));
                    chunk[i].Y = 0;
                }
                splittedContent.Add(chunk);
            }
            return splittedContent;
        }

    }

    public struct FrequencySpectrum
    {

        public readonly Complex[] frequencyDomain;

        private readonly uint samplingFrequency;
        private readonly uint fftWindowSize;


        public FrequencySpectrum(Complex[] frequencyDomain, uint samplingFrequency, uint fftWindowSize)
        {
            if (frequencyDomain.Length == 0)
            {
                throw new ArgumentException("Argument value must be greater than 0", nameof(frequencyDomain));
            }
            if (samplingFrequency == 0)
            {
                throw new ArgumentException("Argument value must be greater than 0", nameof(samplingFrequency));
            }
            if (fftWindowSize == 0)
            {
                throw new ArgumentException("Argument value must be greater than 0", nameof(fftWindowSize));
            }

            this.frequencyDomain = frequencyDomain;
            this.samplingFrequency = samplingFrequency;
            this.fftWindowSize = fftWindowSize;
        }


        //returns magnitude for freq
        public float this[float freq]
        {
            get {
                if (freq >= this.samplingFrequency)
                {
                    throw new IndexOutOfRangeException();
                }

                //find corresponding bin
                float k = freq / ((float)this.samplingFrequency / this.fftWindowSize);
                Complex c = this.frequencyDomain[checked((uint)k)];
                return (float)Math.Sqrt(c.X * c.X + c.Y * c.Y);
            }
        }
    }
}
