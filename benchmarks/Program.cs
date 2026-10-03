using System;
using BenchmarkDotNet.Running;

namespace DotNetMatrix.Benchmarks;

internal static class Program
{
    private static void Main(string[] args)
    {
        if (args.Length == 1 && args[0] == "--validate")
        {
            MultiplicationCandidates.Validate();
            Console.WriteLine("All multiplication candidates satisfy numerical, ownership, and reuse checks.");
            return;
        }
        BenchmarkSwitcher.FromAssembly(typeof(Program).Assembly).Run(args);
    }
}
