using System;
using BenchmarkDotNet.Running;

namespace DotNetMatrix.Benchmarks;

internal static class Program
{
    private static void Main(string[] args)
    {
        if (args.Length == 1 && string.Equals(args[0], "--validate", StringComparison.Ordinal))
        {
            MultiplicationCandidates.Validate();
            Console.WriteLine("All multiplication candidates satisfy numerical, ownership, and reuse checks.");
            return;
        }
        BenchmarkSwitcher.FromAssembly(typeof(Program).Assembly).Run(args);
    }
}
