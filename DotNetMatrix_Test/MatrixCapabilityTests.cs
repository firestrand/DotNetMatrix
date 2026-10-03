using System;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class MatrixCapabilityTests
{
    private static GeneralMatrix M(int rows, params double[] values) => new(rows, values);
    private static void Close(GeneralMatrix expected, GeneralMatrix actual, double tolerance = 1e-11)
    {
        Assert.AreEqual(expected.RowDimension, actual.RowDimension);
        Assert.AreEqual(expected.ColumnDimension, actual.ColumnDimension);
        Assert.IsTrue(expected.Subtract(actual).NormF() <= tolerance * Math.Max(1, expected.NormF()));
    }

    [TestMethod]
    public void PseudoinverseSatisfiesFourIdentitiesAcrossShapesAndDeficientRanks()
    {
        foreach (var a in new[] { M(2, 4, 1, 2, 3), M(3, 1, 2, 3, 4, 5, 7), M(2, 1, 2, 3, 2, 4, 6), new GeneralMatrix(2, 3) })
        {
            var before = a.Copy();
            var p = a.PseudoInverse();
            Close(a, a.Multiply(p).Multiply(a));
            Close(p, p.Multiply(a).Multiply(p));
            Close(a.Multiply(p), a.Multiply(p).Transpose());
            Close(p.Multiply(a), p.Multiply(a).Transpose());
            Close(before, a);
        }
    }

    [TestMethod]
    public void CutoffTiesTruncateAndRankSharesTheSamePolicy()
    {
        var a = GeneralMatrix.Diagonal(1, 5e-5, 1e-4, 2e-4);
        var p = a.PseudoInverse(1e-4);
        Close(GeneralMatrix.Diagonal(1, 0, 0, 5000), p);
        var result = a.SolveMinimumNormWithDiagnostics(GeneralMatrix.Identity(4, 4), 1e-4, true);
        Assert.AreEqual(2, result.NumericalRank);
        Assert.AreEqual(Math.Sqrt(2), result.ResidualNorm, 1e-12);
        Assert.AreEqual(1e-4, a.Subtract(a.Multiply(p).Multiply(a)).Norm2(), 1e-12);
        Assert.AreEqual(20000, result.ConditionNumber!.Value, 1e-6);
        Close(p, result.Solution);
        var truncated = GeneralMatrix.Diagonal(1, 0, 0, 2e-4);
        Close(truncated, a.Multiply(p).Multiply(a));
        Close(truncated, truncated.Multiply(p).Multiply(truncated));
        Close(p, p.Multiply(truncated).Multiply(p));
        Close(truncated.Multiply(p), truncated.Multiply(p).Transpose());
        Close(p.Multiply(truncated), p.Multiply(truncated).Transpose());
        var example = GeneralMatrix.Diagonal(1, 1e-4);
        Close(GeneralMatrix.Diagonal(1, 0), example.PseudoInverse(1e-3));
        Assert.AreEqual(1e-4, example.Subtract(example.Multiply(example.PseudoInverse(1e-3)).Multiply(example)).Norm2(), 1e-12);
        Close(a.PseudoInverse(1e-3), a.Multiply(1000).PseudoInverse(1e-3).Multiply(1000));
        foreach (double cutoff in new[] { -1.0, double.NaN, double.PositiveInfinity })
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => a.PseudoInverse(cutoff));
    }

    [TestMethod]
    public void MinimumNormSolverIsOrthogonalToNullspaceAndHandlesInconsistency()
    {
        var a = M(1, 1, 1);
        var b = M(1, 2, 4);
        var x = a.SolveMinimumNorm(b);
        Close(M(2, 1, 2, 1, 2), x);
        Close(b, a.Multiply(x));
        Close(new GeneralMatrix(1, 2), M(1, 1, -1).Multiply(x));
        var tall = M(3, 1, 0, 0, 1, 0, 0);
        var right = M(3, 2, 3, 4);
        var result = tall.SolveMinimumNormWithDiagnostics(right);
        Close(M(2, 2, 3), result.Solution);
        Assert.AreEqual(4, result.ResidualNorm, 1e-12);
        Assert.AreEqual(2, result.NumericalRank);
        Assert.IsNull(result.ConditionNumber);
        Assert.IsTrue(result.Converged);
        var zero = new GeneralMatrix(2, 3).SolveMinimumNormWithDiagnostics(M(2, 2, 4));
        Assert.AreEqual(0, zero.NumericalRank);
        Assert.AreEqual(Math.Sqrt(20), zero.ResidualNorm, 1e-12);
        Close(new GeneralMatrix(3, 1), zero.Solution);
        Assert.ThrowsExactly<ArgumentException>(() => a.SolveMinimumNorm(new GeneralMatrix(2, 1)));
        Assert.ThrowsExactly<ArgumentNullException>(() => a.SolveMinimumNorm(null!));
    }

    [TestMethod]
    public void MinimumNormAvoidsOverflowFromAnUnrepresentablePseudoinverse()
    {
        foreach (double scale in new[] { 1e-309, 1e200 })
        {
            var a = M(1, scale);
            var result = a.SolveMinimumNormWithDiagnostics(M(1, scale, 2 * scale), 0);
            Close(M(1, 1, 2), result.Solution);
            Assert.AreEqual(0, result.ResidualNorm);
            Assert.AreEqual(1, result.NumericalRank);
        }
        var tall = M(2, 1, 1);
        var large = tall.SolveMinimumNormWithDiagnostics(M(2, double.MaxValue, double.MaxValue), 0);
        Assert.IsTrue(double.IsFinite(large.Solution[0, 0]));
        Assert.AreEqual(1, large.Solution[0, 0] / double.MaxValue, 1e-12);
        Assert.IsTrue(large.ResidualNorm / double.MaxValue < 1e-12);
    }

    [TestMethod]
    public void OrdinaryDiagnosticsAndTrySolveKeepFailureContractsExplicit()
    {
        var a = M(2, 2, 1, 0, 3);
        var b = M(2, 8, 15);
        var result = a.SolveWithDiagnostics(b);
        Close(M(2, 1.5, 5), result.Solution);
        Assert.AreEqual(0, result.ResidualNorm, 1e-12);
        Assert.AreEqual(2, result.NumericalRank);
        Assert.IsNull(result.ConditionNumber);
        Assert.IsTrue(a.SolveWithDiagnostics(b, true).ConditionNumber > 1);
        Assert.IsTrue(a.TrySolve(b, out var x));
        Close(result.Solution, x!);
        Assert.IsFalse(M(2, 1, 2, 2, 4).TrySolve(b, out var missing));
        Assert.IsNull(missing);
        var tall = M(3, 1, 0, 0, 1, 1, 1);
        Assert.IsTrue(tall.TrySolve(M(3, 1, 2, 3), out x));
        Close(M(2, 1, 2), x!);
        Assert.IsFalse(new GeneralMatrix(3, 2).TrySolve(new GeneralMatrix(3, 1), out missing));
        Assert.ThrowsExactly<ArgumentException>(() => a.TrySolve(new GeneralMatrix(1, 1), out _));
        Assert.ThrowsExactly<ArgumentNullException>(() => a.TrySolve(null!, out _));
    }

    [TestMethod]
    public void RepeatedFactorSolvesPreserveInputsAndFactorsWhileBorrowedMutationMatters()
    {
        var a = M(2, 4, 2, 2, 10);
        var aCopy = a.Copy();
        var chol = a.Chol();
        var lu = a.Lud();
        var qr = a.Qrd();
        var tall = M(3, 1, 0, 0, 1, 1, 1);
        var tallQr = tall.Qrd();
        (GeneralMatrix Coefficients, Func<GeneralMatrix, GeneralMatrix> Solve, Func<GeneralMatrix[]> Factors)[] solvers =
        [
            (a, chol.Solve, () => [chol.GetL()]), (a, lu.Solve, () => [lu.L, lu.U]),
            (a, qr.Solve, () => [qr.Q, qr.R, qr.H]), (tall, tallQr.Solve, () => [tallQr.Q, tallQr.R, tallQr.H])
        ];
        var pivot = lu.DoublePivot;
        foreach (var solver in solvers)
        {
            var snapshots = System.Array.ConvertAll(solver.Factors(), factor => factor.Copy());
            foreach (var expected in new[] { M(2, 1, 3), M(2, -2, 1, 4, 2) })
            {
                var b = solver.Coefficients.Multiply(expected);
                var copy = b.Copy();
                Close(expected, solver.Solve(b));
                Close(copy, b, 0);
                var currentFactors = solver.Factors();
                for (int i = 0; i < snapshots.Length; i++) Close(snapshots[i], currentFactors[i], 0);
            }
        }
        CollectionAssert.AreEqual(pivot, lu.DoublePivot);
        Close(aCopy, a, 0);
        var rhs = M(2, 10, 32);
        chol.GetLCopy()[0, 0] = 100;
        Close(M(2, 1, 3), chol.Solve(rhs));
        chol.GetL()[0, 0] = 3;
        Assert.IsFalse(chol.Solve(rhs).ApproximatelyEquals(M(2, 1, 3)));
        Assert.IsTrue(a.Multiply(chol.Solve(rhs)).Subtract(rhs).NormF() > 1);
        Close(M(2, 10, 32), rhs, 0);
    }

    [TestMethod]
    public void ConvenienceMethodsCopyValidateAndPreserveEmptyDimensions()
    {
        var a = M(2, 1, 2, 3, 4, 5, 6);
        Assert.AreEqual(2, a[0, 1]);
        a[0, 1] = 8;
        var row = a.GetRow(0);
        var column = a.GetColumn(1);
        Close(M(1, 1, 8, 3), row);
        Close(M(2, 8, 5), column);
        row[0, 0] = -1;
        column[0, 0] = -1;
        Assert.AreEqual(1, a[0, 0]);
        Assert.AreEqual(8, a[0, 1]);
        Close(M(3, 2, 0, 0, 0, 3, 0, 0, 0, 4), GeneralMatrix.Diagonal(2, 3, 4));
        Close(M(2, 1, 8, 3, 1, 8, 3), row = M(1, 1, 8, 3).ConcatenateRows(M(1, 1, 8, 3)));
        Close(a, a.GetColumn(0).ConcatenateColumns(a.GetColumn(1)).ConcatenateColumns(a.GetColumn(2)));
        var emptyColumn = new GeneralMatrix(0, 3).GetColumn(1);
        Assert.AreEqual(0, emptyColumn.RowDimension);
        Assert.AreEqual(1, emptyColumn.ColumnDimension);
        Assert.AreEqual(0, new GeneralMatrix(2, 0).GetRow(0).ColumnDimension);
        Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => a.GetRow(-1));
        Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => a.GetColumn(3));
        Assert.ThrowsExactly<ArgumentException>(() => a.ConcatenateRows(new GeneralMatrix(1, 2)));
        Assert.ThrowsExactly<ArgumentException>(() => a.ConcatenateColumns(new GeneralMatrix(1, 1)));
        Assert.ThrowsExactly<ArgumentNullException>(() => a.ConcatenateRows(null!));
        Assert.ThrowsExactly<ArgumentNullException>(() => GeneralMatrix.Diagonal(null!));
        var buffer = new double[8];
        a.CopyTo(buffer.AsSpan(1, 6));
        CollectionAssert.AreEqual(new[] { 0.0, 1, 8, 3, 4, 5, 6, 0 }, buffer);
        a.CopyTo(buffer.AsSpan(1, 6), true);
        CollectionAssert.AreEqual(new[] { 0.0, 1, 4, 8, 5, 3, 6, 0 }, buffer);
        Assert.ThrowsExactly<ArgumentException>(() => a.CopyTo(new double[5]));
    }

    [TestMethod]
    public void ApproximateComparisonHasExplicitToleranceAndSpecialValueRules()
    {
        var a = M(1, 1e10, 0, double.PositiveInfinity);
        Assert.IsTrue(a.ApproximatelyEquals(M(1, 1e10 + 1, 1e-10, double.PositiveInfinity), 1e-9, 1e-9));
        Assert.IsFalse(a.ApproximatelyEquals(M(1, 1e10 + 1, 1, double.NegativeInfinity)));
        Assert.IsFalse(M(1, double.NaN).ApproximatelyEquals(M(1, double.NaN)));
        Assert.IsTrue(M(1, -0.0).ApproximatelyEquals(M(1, 0.0), 0, 0));
        Assert.IsFalse(a.ApproximatelyEquals(null));
        Assert.IsFalse(a.ApproximatelyEquals(new GeneralMatrix(3, 1)));
        Assert.IsFalse(M(1, double.PositiveInfinity).ApproximatelyEquals(M(1, 1)));
        foreach (double invalid in new[] { -1.0, double.NaN, double.PositiveInfinity })
        {
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => a.ApproximatelyEquals(a, invalid));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => a.ApproximatelyEquals(a, 0, invalid));
        }
        // Finite opposite extremes must not compare equal through infinity overflow.
        Assert.IsFalse(M(1, double.MaxValue).ApproximatelyEquals(M(1, -double.MaxValue), 0, 1));
    }

    [TestMethod]
    public void CallerRandomIsRepeatableAndDoesNotBorrowStorage()
    {
        var first = GeneralMatrix.Random(2, 3, new Random(42));
        var second = GeneralMatrix.Random(2, 3, new Random(42));
        Close(first, second, 0);
        foreach (double value in first.RowPackedCopy) Assert.IsTrue(value >= 0 && value < 1);
        second[0, 0] = -1;
        Assert.IsTrue(first[0, 0] >= 0);
        Assert.ThrowsExactly<ArgumentNullException>(() => GeneralMatrix.Random(1, 1, null!));
        Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GeneralMatrix.Random(-1, 1, new Random()));
    }
}
