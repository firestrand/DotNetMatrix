using System;
using System.Linq;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class EigenvalueTests
{
    [TestMethod]
    public void SymmetricMatricesHaveSortedRealEigenvaluesAndOrthogonalEigenvectors()
    {
        GeneralMatrix[] matrices =
        [
            new(new[] { new[] { -7.0 } }),
            new(4, 4),
            GeneralMatrix.Identity(5, 5).Multiply(3),
            new(new[] { new[] { 4.0, 0, 0 }, new[] { 0.0, -2, 0 }, new[] { 0.0, 0, 4 } }),
            new(new[] { new[] { 2.0, 1 }, new[] { 1.0, 2 } }),
            new(new[] { new[] { 2.0, -1, 0 }, new[] { -1.0, 2, -1 }, new[] { 0.0, -1, 2 } }),
            new(new[] { new[] { 1.0, 2, 3 }, new[] { 2.0, 4, 6 }, new[] { 3.0, 6, 9 } }),
            new(new[] { new[] { 0.0, -2, 4 }, new[] { -2.0, 0, -3 }, new[] { 4.0, -3, 0 } })
        ];

        foreach (GeneralMatrix matrix in matrices)
        {
            AssertSymmetricDecomposition(matrix);
        }

        var simple = matrices[4].Eigen();
        Assert.AreEqual(1.0, simple.RealEigenvalues[0], 1e-12);
        Assert.AreEqual(3.0, simple.RealEigenvalues[1], 1e-12);
        var rankOne = matrices[6].Eigen();
        Assert.AreEqual(0.0, rankOne.RealEigenvalues[0], 1e-12);
        Assert.AreEqual(0.0, rankOne.RealEigenvalues[1], 1e-12);
        Assert.AreEqual(14.0, rankOne.RealEigenvalues[2], 1e-12);
    }

    [TestMethod]
    public void SymmetricDenseAndScaledMatricesReconstructWithinRelativePrecision()
    {
        var random = new Random(1729);
        for (int size = 2; size <= 14; size++)
        {
            double[][] entries = Allocate(size);
            for (int row = 0; row < size; row++)
            {
                for (int column = 0; column <= row; column++)
                {
                    entries[row][column] = entries[column][row] = random.NextDouble() * 20 - 10;
                }
            }

            var matrix = new GeneralMatrix(entries);
            AssertSymmetricDecomposition(matrix);
            AssertSymmetricDecomposition(matrix.Multiply(1e-100));
            AssertSymmetricDecomposition(matrix.Multiply(1e100));
        }
    }

    [TestMethod]
    public void TriangularMatricesRetainKnownEigenvaluesIncludingRepeatedRoots()
    {
        GeneralMatrix[] matrices =
        [
            new(new[] { new[] { 1.0, 4, -3 }, new[] { 0.0, -2, 7 }, new[] { 0.0, 0, 5 } }),
            new(new[] { new[] { 2.0, 1, 0 }, new[] { 0.0, 2, 1 }, new[] { 0.0, 0, 2 } }),
            new(new[] { new[] { 0.0, 1, 0 }, new[] { 0.0, 0, 1 }, new[] { 0.0, 0, 0 } }),
            new(new[] { new[] { 1.0, 0, 0 }, new[] { 4.0, -2, 0 }, new[] { -3.0, 7, 5 } })
        ];

        foreach (GeneralMatrix matrix in matrices)
        {
            EigenvalueDecomposition eigen = AssertEigensystem(matrix);
            double[] expected = Enumerable.Range(0, matrix.RowDimension)
                .Select(index => matrix.GetElement(index, index)).OrderBy(value => value).ToArray();
            double[] actual = eigen.RealEigenvalues.OrderBy(value => value).ToArray();
            for (int index = 0; index < expected.Length; index++)
            {
                Assert.AreEqual(expected[index], actual[index], 1e-10);
                Assert.AreEqual(0.0, eigen.ImagEigenvalues[index], 1e-12);
            }
        }
    }

    [TestMethod]
    public void RotationBlocksRepresentComplexConjugateEigenvaluesWithCorrectSigns()
    {
        GeneralMatrix[] matrices =
        [
            new(new[] { new[] { 0.0, -1 }, new[] { 1.0, 0 } }),
            new(new[] { new[] { 3.0, 2 }, new[] { -8.0, 3 } }),
            new(new[] { new[] { -2.0, -8 }, new[] { 2.0, -2 } }),
            new(new[] { new[] { 0.0, -1, 4, 2 }, new[] { 1.0, 0, 3, 5 }, new[] { 0.0, 0, 3, -2 }, new[] { 0.0, 0, 2, 3 } }),
            new(new[] { new[] { 0.0, -1, 1, 0 }, new[] { 1.0, 0, 0, 1 }, new[] { 0.0, 0, 0, -1 }, new[] { 0.0, 0, 1, 0 } }),
            new(new[] { new[] { 3.0, -2, 7 }, new[] { 2.0, 3, -4 }, new[] { 0.0, 0, -5 } })
        ];

        foreach (GeneralMatrix matrix in matrices)
        {
            EigenvalueDecomposition eigen = AssertEigensystem(matrix);
            Assert.IsTrue(eigen.ImagEigenvalues.Any(value => value > 0));
            AssertComplexBlocks(eigen);
        }

        EigenvalueDecomposition rotation = matrices[0].Eigen();
        Assert.AreEqual(0.0, rotation.RealEigenvalues[0], 1e-12);
        Assert.AreEqual(0.0, rotation.RealEigenvalues[1], 1e-12);
        Assert.AreEqual(1.0, rotation.ImagEigenvalues[0], 1e-12);
        Assert.AreEqual(-1.0, rotation.ImagEigenvalues[1], 1e-12);
        EigenvalueDecomposition shifted = matrices[1].Eigen();
        Assert.AreEqual(3.0, shifted.RealEigenvalues[0], 1e-12);
        Assert.AreEqual(4.0, shifted.ImagEigenvalues[0], 1e-12);
    }

    [TestMethod]
    public void NonsymmetricDenseMatricesSatisfyEigenvectorAndSpectralIdentities()
    {
        var random = new Random(271828);
        for (int size = 2; size <= 14; size++)
        {
            for (int fixture = 0; fixture < 8; fixture++)
            {
                double[][] entries = Allocate(size);
                for (int row = 0; row < size; row++)
                {
                    for (int column = 0; column < size; column++)
                    {
                        entries[row][column] = random.NextDouble() * 20 - 10;
                    }
                }

                GeneralMatrix matrix = new(entries);
                EigenvalueDecomposition eigen = AssertEigensystem(matrix);
                AssertComplexBlocks(eigen);
                // tr(A^2) = sum(lambda^2), including complex conjugate pairs.
                double spectralSquareTrace = eigen.RealEigenvalues.Select((real, index) =>
                    real * real - eigen.ImagEigenvalues[index] * eigen.ImagEigenvalues[index]).Sum();
                Assert.AreEqual(matrix.Multiply(matrix).Trace(), spectralSquareTrace,
                    1e-9 * Math.Max(1, matrix.Norm1() * matrix.Norm1()));
            }
        }
    }

    [TestMethod]
    public void HessenbergAndNearlyDefectiveMatricesPreserveTheEigenvectorEquation()
    {
        GeneralMatrix[] matrices =
        [
            new(new[] { new[] { 0.0, 1, 0, 0 }, new[] { 1.0, 0, 2e-7, 0 }, new[] { 0.0, -2e-7, 0, 1 }, new[] { 0.0, 0, 1, 0 } }),
            new(new[] { new[] { 0.0, 0, 0, -1 }, new[] { 1.0, 0, 0, 0 }, new[] { 0.0, 1, 0, 0 }, new[] { 0.0, 0, 1, 0 } }),
            new(new[] { new[] { 0.0, 0, -6 }, new[] { 1.0, 0, 11 }, new[] { 0.0, 1, -6 } }),
            new(new[] { new[] { 1.0, 1, 0 }, new[] { 1e-12, 1, 1 }, new[] { 0.0, 1e-12, 1 } }),
            new(new[] { new[] { 1.0, 1e10, 0 }, new[] { 1e-10, 2, 1e10 }, new[] { 0.0, -1e-10, 3 } })
        ];

        foreach (GeneralMatrix matrix in matrices)
        {
            AssertComplexBlocks(AssertEigensystem(matrix));
        }

        foreach (double scale in new[] { 1e-100, 1e100 })
        {
            AssertEigensystem(matrices[1].Multiply(scale));
        }
    }

    [TestMethod]
    public void DecompositionDoesNotMutateTheInputMatrix()
    {
        var matrix = new GeneralMatrix(new[] { new[] { 1.0, 2, 3 }, new[] { -4.0, 5, 6 }, new[] { 7.0, 8, -9 } });
        double[] original = matrix.RowPackedCopy;
        _ = matrix.Eigen();
        CollectionAssert.AreEqual(original, matrix.RowPackedCopy);
    }

    private static void AssertSymmetricDecomposition(GeneralMatrix matrix)
    {
        EigenvalueDecomposition eigen = AssertEigensystem(matrix);
        GeneralMatrix vectors = eigen.GetV();
        Assert.IsTrue(vectors.Transpose().Multiply(vectors)
            .Subtract(GeneralMatrix.Identity(matrix.RowDimension, matrix.RowDimension)).Norm1() < 1e-10);
        double relativeError = vectors.Multiply(eigen.D).Multiply(vectors.Transpose())
            .Subtract(matrix).Norm1() / Math.Max(matrix.Norm1(), 1e-300);
        Assert.IsTrue(relativeError < 1e-10, $"Symmetric reconstruction residual: {relativeError:R}");
        for (int index = 0; index < eigen.RealEigenvalues.Length; index++)
        {
            Assert.AreEqual(0.0, eigen.ImagEigenvalues[index]);
            if (index > 0)
            {
                Assert.IsTrue(eigen.RealEigenvalues[index - 1] <= eigen.RealEigenvalues[index]);
            }
        }
    }

    private static EigenvalueDecomposition AssertEigensystem(GeneralMatrix matrix)
    {
        EigenvalueDecomposition eigen = matrix.Eigen();
        GeneralMatrix vectors = eigen.GetV();
        Assert.AreEqual(matrix.RowDimension, vectors.RowDimension);
        Assert.AreEqual(matrix.ColumnDimension, vectors.ColumnDimension);
        Assert.IsTrue(eigen.RealEigenvalues.All(double.IsFinite));
        Assert.IsTrue(eigen.ImagEigenvalues.All(double.IsFinite));
        Assert.IsTrue(vectors.RowPackedCopy.All(double.IsFinite));
        Assert.IsTrue(vectors.Norm1() > 0);
        double denominator = Math.Max(matrix.Norm1() * vectors.Norm1(), 1e-300);
        double residual = matrix.Multiply(vectors).Subtract(vectors.Multiply(eigen.D)).Norm1() / denominator;
        Assert.IsTrue(residual < 1e-10, $"A*V = V*D relative residual: {residual:R}");
        Assert.AreEqual(matrix.Trace(), eigen.RealEigenvalues.Sum(),
            Math.Max(matrix.Norm1() * 1e-10, 1e-300));
        return eigen;
    }

    private static void AssertComplexBlocks(EigenvalueDecomposition eigen)
    {
        GeneralMatrix diagonal = eigen.D;
        for (int index = 0; index < eigen.RealEigenvalues.Length; index++)
        {
            Assert.AreEqual(eigen.RealEigenvalues[index], diagonal.GetElement(index, index));
            double imaginary = eigen.ImagEigenvalues[index];
            if (imaginary > 0)
            {
                Assert.IsTrue(index + 1 < eigen.RealEigenvalues.Length);
                Assert.AreEqual(imaginary, -eigen.ImagEigenvalues[index + 1]);
                Assert.AreEqual(eigen.RealEigenvalues[index], eigen.RealEigenvalues[index + 1]);
                Assert.AreEqual(imaginary, diagonal.GetElement(index, index + 1));
            }
            else if (imaginary < 0)
            {
                Assert.IsTrue(index > 0);
                Assert.AreEqual(imaginary, diagonal.GetElement(index, index - 1));
            }
        }
    }

    private static double[][] Allocate(int size) =>
        Enumerable.Range(0, size).Select(_ => new double[size]).ToArray();
}
