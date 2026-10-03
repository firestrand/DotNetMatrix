using System;

namespace DotNetMatrix;

/// <summary>A numerical decomposition exhausted its configured iteration budget.</summary>
public sealed class DecompositionConvergenceException : SystemException
{
    /// <summary>Creates an error identifying the algorithm and exhausted budget.</summary>
    public DecompositionConvergenceException(string algorithm, int maxIterations)
        : base($"{algorithm} did not converge within {maxIterations} iterations.")
    {
        Algorithm = algorithm;
        MaxIterations = maxIterations;
    }

    /// <summary>The algorithm that failed to converge.</summary>
    public string Algorithm { get; }

    /// <summary>The configured total iteration budget.</summary>
    public int MaxIterations { get; }
}
