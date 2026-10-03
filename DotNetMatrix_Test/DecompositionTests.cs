using System;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class DecompositionTests
{
    [TestMethod]
    public void CholeskyReconstructsPositiveDefiniteMatrixWithNonunitDiagonalFactor()
    {
        var lower = Matrix([2, 0, 0], [-1, 3, 0], [4, 2, 5]);
        var input = lower.Multiply(lower.Transpose());
        var decomposition = input.Chol();
        Assert.IsTrue(decomposition.SPD);
        AssertMatrix(lower, decomposition.GetL());
        AssertMatrix(input, decomposition.GetL().Multiply(decomposition.GetL().Transpose()));
    }

    [TestMethod]
    public void CholeskySolvesMultipleRightHandSidesWithUnitDiagonalFactor()
    {
        var lower = Matrix([1, 0, 0], [-1, 1, 0], [4, 2, 1]);
        var input = lower.Multiply(lower.Transpose());
        var decomposition = input.Chol();
        var expected = Matrix([1, -2], [3, 4], [-1, 5]);
        var right = input.Multiply(expected);
        var original = right.Copy();
        AssertMatrix(expected, decomposition.Solve(right));
        AssertMatrix(original, right);
    }

    [TestMethod]
    public void CholeskyNonunitDiagonalSolveSatisfiesEquation()
    {
        var input = Matrix([4, 2], [2, 10]);
        var right = Matrix([10, 0], [32, 36]);
        var original = right.Copy();
        var solution = input.Chol().Solve(right);
        AssertMatrix(Matrix([1, -2], [3, 4]), solution);
        AssertMatrix(right, input.Multiply(solution));
        AssertMatrix(original, right);
    }

    [TestMethod]
    public void CholeskyRejectsNonSymmetricIndefiniteSingularAndRectangularMatrices()
    {
        GeneralMatrix[] invalid =
        [
            Matrix([2, 1], [0, 2]),
            Matrix([1, 2], [2, 1]),
            Matrix([0, 0], [0, 1])
        ];
        foreach (var input in invalid)
        {
            var decomposition = new CholeskyDecomposition(input);
            Assert.IsFalse(decomposition.SPD);
            Assert.ThrowsExactly<SystemException>(() => decomposition.Solve(GeneralMatrix.Identity(2, 2)));
        }
        Assert.ThrowsExactly<ArgumentException>(() => new CholeskyDecomposition(Matrix([1, 0, 5], [0, 1, 6])));
        Assert.ThrowsExactly<ArgumentException>(() => Matrix([1]).Chol().Solve(new GeneralMatrix(2, 1)));
    }

    [TestMethod]
    public void LuPivotingReconstructsMatrixAndSolvesWithoutMutatingInputs()
    {
        var input = Matrix([0, 2, 1], [3, -1, 2], [1, 1, 4]);
        var original = input.Copy();
        var decomposition = new LUDecomposition(input);
        Assert.IsTrue(decomposition.IsNonSingular);
        AssertMatrix(input.GetMatrix(decomposition.Pivot, 0, 2), decomposition.L.Multiply(decomposition.U));
        Assert.AreEqual(-16d, decomposition.Determinant(), 1e-12);
        var pivots = decomposition.Pivot;
        var doublePivots = decomposition.DoublePivot;
        for (var i = 0; i < pivots.Length; i++)
        {
            Assert.AreEqual((double)pivots[i], doublePivots[i]);
        }
        pivots[0] = -1;
        Assert.AreNotEqual(-1, decomposition.Pivot[0]);
        var expected = Matrix([1, -2], [3, 4], [-1, 5]);
        var right = input.Multiply(expected);
        var rightCopy = right.Copy();
        AssertMatrix(expected, decomposition.Solve(right));
        AssertMatrix(original, input);
        AssertMatrix(rightCopy, right);
    }

    [TestMethod]
    public void LuTallAndSingularMatricesRetainFactorizationAndErrorContracts()
    {
        var tall = Matrix([3, 0], [1, 2], [-1, 4], [0, 1]);
        var decomposition = new LUDecomposition(tall);
        AssertMatrix(tall.GetMatrix(decomposition.Pivot, 0, 1), decomposition.L.Multiply(decomposition.U));
        Assert.ThrowsExactly<ArgumentException>(() => decomposition.Determinant());
        Assert.ThrowsExactly<ArgumentException>(() => decomposition.Solve(new GeneralMatrix(3, 1)));
        var singular = new LUDecomposition(Matrix([1, 2], [2, 4]));
        Assert.IsFalse(singular.IsNonSingular);
        Assert.AreEqual(0d, singular.Determinant());
        Assert.ThrowsExactly<SystemException>(() => singular.Solve(new GeneralMatrix(2, 1)));
        var diagonal = new LUDecomposition(Matrix([2, 0], [0, -3]));
        Assert.AreEqual(-6d, diagonal.Determinant());
        AssertMatrix(GeneralMatrix.Identity(2, 2), diagonal.L);
    }

    [TestMethod]
    public void QrReconstructsTallMatrixWithOrthogonalColumnsAndHouseholderVectors()
    {
        var input = Matrix([-2, 3, 1], [4, -1, 2], [1, 5, -2], [3, 0, 4], [0, 2, 1]);
        var original = input.Copy();
        var decomposition = new QRDecomposition(input);
        Assert.IsTrue(decomposition.FullRank);
        AssertMatrix(input, decomposition.Q.Multiply(decomposition.R));
        AssertMatrix(GeneralMatrix.Identity(3, 3), decomposition.Q.Transpose().Multiply(decomposition.Q));
        var h = decomposition.H;
        for (var i = 0; i < h.RowDimension; i++)
        {
            for (var j = 0; j < h.ColumnDimension; j++)
            {
                if (i < j)
                {
                    Assert.AreEqual(0d, h.GetElement(i, j));
                }
            }
        }
        AssertMatrix(original, input);
    }

    [TestMethod]
    public void QrLeastSquaresHasKnownSolutionAndOrthogonalResidual()
    {
        var input = Matrix([1, 0], [1, 1], [1, 2], [1, 3]);
        var right = Matrix([1, 0], [2, 1], [2, 4], [4, 9]);
        var original = right.Copy();
        var solution = new QRDecomposition(input).Solve(right);
        AssertMatrix(Matrix([0.9, -1], [0.9, 3]), solution);
        var residual = input.Multiply(solution).Subtract(right);
        AssertMatrix(new GeneralMatrix(2, 2), input.Transpose().Multiply(residual));
        AssertMatrix(original, right);
    }

    [TestMethod]
    public void QrRankDeficiencyAllowsReconstructionButRejectsSolve()
    {
        var input = Matrix([1, 0], [0, 0], [0, 0]);
        var decomposition = new QRDecomposition(input);
        Assert.IsFalse(decomposition.FullRank);
        AssertMatrix(input, decomposition.Q.Multiply(decomposition.R));
        Assert.ThrowsExactly<SystemException>(() => decomposition.Solve(new GeneralMatrix(3, 1)));
        Assert.ThrowsExactly<ArgumentException>(() => decomposition.Solve(new GeneralMatrix(2, 1)));
        var zero = new QRDecomposition(new GeneralMatrix(3, 2));
        Assert.IsFalse(zero.FullRank);
        AssertMatrix(new GeneralMatrix(3, 2), zero.H);
        AssertMatrix(new GeneralMatrix(3, 2), zero.Q.Multiply(zero.R));
    }

    [TestMethod]
    public void SvdReconstructsDiverseDeterministicMatricesWithOrthogonalVectors()
    {
        GeneralMatrix[] fixtures =
        [
            Matrix([1, 2, 3], [4, 5, 6]),
            Matrix([1, 1]),
            new GeneralMatrix(2, 5),
            Matrix([3]),
            Matrix([-2]),
            Matrix([0]),
            Matrix([1, 0, 0], [0, 5, 0], [0, 0, -3]),
            Matrix([0, 2, 0], [3, 0, 0], [0, 0, 1]),
            Matrix([1, 2, 3], [2, 4, 6], [0, 0, 0]),
            Matrix([0, 1, -2], [0, 3, 4], [0, -2, 1], [0, 1, 3]),
            new GeneralMatrix(4, 3)
        ];
        foreach (var fixture in fixtures)
        {
            AssertSvd(fixture);
        }
        // Dense fixtures exercise iterative convergence rather than only diagonal cases.
        for (var columns = 2; columns <= 8; columns++)
        {
            var fixture = new GeneralMatrix(columns + 2, columns);
            for (var i = 0; i < fixture.RowDimension; i++)
            {
                for (var j = 0; j < columns; j++)
                {
                    fixture.SetElement(i, j, Math.Sin((i + 1) * (j + 2)) + (i == j ? 2 : 0));
                }
            }
            AssertSvd(fixture);
        }
    }

    [TestMethod]
    public void SvdReportsKnownSingularValuesConditionNumberAndNumericalRank()
    {
        var input = Matrix([3, 0, 0], [0, -5, 0], [0, 0, 1]);
        var decomposition = new SingularValueDecomposition(input);
        CollectionAssert.AreEqual(new double[] { 5, 3, 1 }, decomposition.SingularValues);
        Assert.AreEqual(5d, decomposition.Norm2());
        Assert.AreEqual(5d, decomposition.Condition());
        Assert.AreEqual(3, decomposition.Rank());
        Assert.AreEqual(5d, input.Norm2());
        Assert.AreEqual(5d, input.Condition());
        Assert.AreEqual(3, input.Rank());
        var singular = new SingularValueDecomposition(Matrix([2, 0], [0, 0]));
        Assert.AreEqual(1, singular.Rank());
        Assert.IsTrue(double.IsPositiveInfinity(singular.Condition()));
        var tiny = new SingularValueDecomposition(Matrix([1, 0], [0, 1e-18]));
        Assert.AreEqual(1, tiny.Rank());
        var zero = new SingularValueDecomposition(new GeneralMatrix(2, 2));
        Assert.AreEqual(0, zero.Rank());
        Assert.AreEqual(0d, zero.Norm2());
        Assert.IsTrue(double.IsNaN(zero.Condition()));
    }

    [TestMethod]
    public void FrobeniusNormAvoidsOverflowAndUnderflowAndHandlesZerosAndSigns()
    {
        Assert.AreEqual(5d, Matrix([-3, 4]).NormF(), 1e-12);
        Assert.AreEqual(5d, Matrix([4, -3]).NormF(), 1e-12);
        Assert.AreEqual(0d, new GeneralMatrix(2, 2).NormF());
        Assert.AreEqual(5d, Matrix([0, -5, 0]).NormF(), 1e-12);
        Assert.AreEqual(5d, Matrix([3e200, 4e200]).NormF() / 1e200, 1e-12);
        Assert.AreEqual(5d, Matrix([3e-200, 4e-200]).NormF() / 1e-200, 1e-12);
    }

    private static GeneralMatrix Matrix(params double[][] rows) => new(rows);

    private static void AssertSvd(GeneralMatrix input)
    {
        var original = input.Copy();
        var decomposition = new SingularValueDecomposition(input);
        var u = decomposition.GetU();
        var v = decomposition.GetV();
        AssertMatrix(input, u.Multiply(decomposition.S).Multiply(v.Transpose()));
        AssertMatrix(GeneralMatrix.Identity(Math.Min(input.RowDimension, input.ColumnDimension), Math.Min(input.RowDimension, input.ColumnDimension)), u.Transpose().Multiply(u));
        AssertMatrix(GeneralMatrix.Identity(Math.Min(input.RowDimension, input.ColumnDimension), Math.Min(input.RowDimension, input.ColumnDimension)), v.Transpose().Multiply(v));
        var values = decomposition.SingularValues;
        for (var i = 0; i < values.Length; i++)
        {
            Assert.IsTrue(values[i] >= 0);
            if (i > 0)
            {
                Assert.IsTrue(values[i - 1] >= values[i]);
            }
        }
        AssertMatrix(original, input);
    }

    private static void AssertMatrix(GeneralMatrix expected, GeneralMatrix actual)
    {
        Assert.AreEqual(expected.RowDimension, actual.RowDimension);
        Assert.AreEqual(expected.ColumnDimension, actual.ColumnDimension);
        for (var i = 0; i < expected.RowDimension; i++)
        {
            for (var j = 0; j < expected.ColumnDimension; j++)
            {
                var value = expected.GetElement(i, j);
                Assert.AreEqual(value, actual.GetElement(i, j), 1e-10 * Math.Max(1, Math.Abs(value)), $"Element [{i},{j}]");
            }
        }
    }
}
