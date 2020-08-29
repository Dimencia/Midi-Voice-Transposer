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

namespace Midi_Transposer_Console
{
    /// <summary>
    /// Handles conversion of sound files to MIDI, each instance of a Converter has its own settings etc
    /// </summary>
    public class Converter
    {
        private readonly float[] MidiFrequencyMap;

        public Converter()
        {
            MidiFrequencyMap = GetMidiFrequencyMap();
        }

        public void ConvertWavToMidi(string inputFile)
        {
            // 4096 is about 40ms, 2048 is about 20ms, good for speech
            uint windowSize = 4096;
            float noteThreshold = 0.1f; // Amplitudes above this get translated to a midi note
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
            int maxFreqToGraph = 12000;
            float xAxisRatio = (float)imageSize / timesteps;
            float yAxisRatio = (float)imageSize / maxFreqToGraph;
            float maxAmp = 0;

            using (Graphics g = Graphics.FromImage(b))
            {
                g.FillRectangle(Brushes.Black, 0, 0, imageSize, imageSize);
                for (int i = 0; i < timesteps; i++)
                {
                    for (uint j = 0; j < maxFreqToGraph; j++)
                    {
                        var amplitude = spectrumList[i][j];
                        if (amplitude > maxAmp)
                            maxAmp = amplitude;
                        amplitude *= 25; // Scale it up so you can see anything for basic display
                        // While running we're calculating the amp normalization for the MIDI to use
                        float x = i * xAxisRatio;
                        float y = imageSize - (j * yAxisRatio);
                        //if (amplitude > noteThreshold)
                        //    Console.WriteLine(amplitude);
                        int color = (int)(amplitude * 255);
                        if (color > 255)
                            color = 255;
                        if (color < 0)
                            color = 0;

                        // TODO: Change this all to lockbits stuff if it ever becomes a perf issue.  Or just remove it.
                        if (color > 0)
                            g.DrawLine(new Pen(Color.FromArgb(color, color, color), 1), x, y, x + 1, y + 1);
                    }
                }
                g.Save();
            }
            float ampMult = 1.0f / maxAmp;
            Console.WriteLine($"Amp Mult: {ampMult}");

            b.Save("test.png", ImageFormat.Png);

            Console.WriteLine("PNG generated, creating midi");
            var midiFile = new MidiFile();
            TempoMap tempoMap = midiFile.GetTempoMap();

            var trackChunk = new TrackChunk();
            using (var notesManager = trackChunk.ManageNotes())
            {
                NotesCollection notes = notesManager.Notes;
                for (int i = 0; i < timesteps; i++)
                {
                    foreach (uint j in MidiFrequencyMap)
                    {
                        var amplitude = Math.Clamp(spectrumList[i][j] * ampMult * 4,0,1);
                        // I tried taking the log of amp (and even 10^amp), but it just kept spamming 1's in either case.
                        // But really, the amps are so low, something seems wrong.
                        // ampMult is often nearly 10.  And often, somehow, no values end up anywhere near 1.0 amplitude even after applying it
                        // And basically no matter what I do, it's not loud enough


                        // Something's not right here, other than that, I think.
                        // It seems high notes get a lot more volume than low ones
                        // Which sounds like a log scale human hearing thing
                        // ... Or just the piano soundfont being weird
                        if (amplitude > noteThreshold)
                        {
                            Console.WriteLine(amplitude);
                            // Convert i into an absolute time...
                            // Which means, taking the length of the wav/i
                            
                            double durationMs = (wavDuration / timesteps).TotalMilliseconds;
                            double positionMs = i * durationMs;
                            var note = new Melanchall.DryWetMidi.Interaction.Note(
                        SevenBitNumber.Parse(GetMidiNoteNumber(j).ToString()),
                        LengthConverter.ConvertFrom(
                            new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)durationMs),
                            0,
                            tempoMap),
                        LengthConverter.ConvertFrom(
                            new MetricTimeSpan(hours: 0, minutes: 0, seconds: 0, milliseconds: (int)positionMs),
                            0,
                            tempoMap));
                            note.Velocity = (SevenBitNumber)(amplitude * 127);
                            notes.Add(note);
                        }
                    }
                }
            }

            midiFile.Chunks.Add(trackChunk);
            midiFile.Write("Vocal_Conversion.mid", true);
            // TODO: Use Polyphone (standalone app) to create a soundfont that is a plain tone
            // And either package this with the midi somehow, or just, use it for playback in midieditor or wherever
        }


        private int GetMidiNoteNumber(uint frequency)
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
            if (mid == a.Length-1 || Math.Abs(item - a[mid + 1]) > Math.Abs(item - a[mid]))
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
            for(int i = 0; i < 128; i++)
            {
                result[i] = 8.1758f * (float)(Math.Pow(2, i/12f)); // Internet said this is how to get them
            }
            return result;
        }

        private FrequencySpectrum ToFrequencySpectrum(Complex[] timeDomain, uint sampleRate, uint windowSize)
        {
            int m = (int)Math.Log(timeDomain.Length, 2); // So, either timeDomain.Length needs to be a power of 2
            // Or this needs to be changed to take windowSize instead.  Not sure

            //true = forward fft
            FastFourierTransform.FFT(true, m, timeDomain);
            return new FrequencySpectrum(timeDomain, sampleRate, windowSize);
        }

        public IList<FrequencySpectrum> Fft(uint windowSize, float[] audioContent, uint sampleRate)
        {
            IList<Complex[]> timeDomainChunks = this.SplitInChunks(audioContent, windowSize);
            return timeDomainChunks.Select(x => ToFrequencySpectrum(x, sampleRate, windowSize)).ToList();
        }

        private IList<Complex[]> SplitInChunks(float[] audioContent, uint chunkSize)
        {
            IList<Complex[]> splittedContent = new List<Complex[]>();

            for (uint k = 0; k < audioContent.Length; k += chunkSize)
            {
                long size = k + chunkSize < audioContent.Length ? chunkSize : audioContent.Length - k;
                Complex[] chunk = new Complex[size];

                for (int i = 0; i < chunk.Length; i++)
                {
                    //i've tried windowing here but didn't seem to help me
                    chunk[i].X = audioContent[k + i];
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
        public float this[uint freq]
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
