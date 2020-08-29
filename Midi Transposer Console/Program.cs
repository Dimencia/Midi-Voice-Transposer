using System;

namespace Midi_Transposer_Console
{
    class Program
    {
        static void Main(string[] args)
        {
            new Converter().ConvertWavToMidi(@"C:\Users\Dimen_hmzu9w6\Downloads\bohemian_vocal_snippet.wav");
        }
    }
}
