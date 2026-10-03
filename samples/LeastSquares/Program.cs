using System;
using DotNetMatrix;

// An analytic example: y = 2 + 3*x, not an asserted real-world dataset.
var design = new GeneralMatrix(4, new[] { 1.0, 0, 1, 1, 1, 2, 1, 3 });
var observations = new GeneralMatrix(4, new[] { 2.0, 5, 8, 11 });
var fit = design.SolveMinimumNormWithDiagnostics(observations);
if (fit.NumericalRank != 2 || fit.ResidualNorm > 1e-10 ||
    Math.Abs(fit.Solution[0, 0] - 2) > 1e-10 || Math.Abs(fit.Solution[1, 0] - 3) > 1e-10)
    throw new InvalidOperationException("Known coefficient fit did not satisfy its analytic oracle.");
Console.WriteLine($"Intercept: {fit.Solution[0, 0]:G6}; slope: {fit.Solution[1, 0]:G6}; rank: {fit.NumericalRank}; residual: {fit.ResidualNorm:G6}");

// Duplicate columns: many coefficients fit; the minimum-norm choice splits them evenly.
var deficient = new GeneralMatrix(3, new[] { 1.0, 1, 1, 1, 1, 1 });
var deficientFit = deficient.SolveMinimumNormWithDiagnostics(new GeneralMatrix(3, new[] { 4.0, 4, 4 }));
if (deficientFit.NumericalRank != 1 || deficientFit.ResidualNorm > 1e-10 ||
    Math.Abs(deficientFit.Solution[0, 0] - 2) > 1e-10 || Math.Abs(deficientFit.Solution[1, 0] - 2) > 1e-10)
    throw new InvalidOperationException("Rank-deficient fit did not satisfy minimum-norm oracle.");
Console.WriteLine($"Duplicate-column coefficients: {deficientFit.Solution[0, 0]:G6}, {deficientFit.Solution[1, 0]:G6}; rank: {deficientFit.NumericalRank}");
