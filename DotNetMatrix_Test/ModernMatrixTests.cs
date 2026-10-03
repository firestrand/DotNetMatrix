using System;
using System.Runtime.Serialization;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class ModernMatrixTests
{
    private static GeneralMatrix Sample() => new(new[] { new[] { 2.0, -4.0 }, new[] { 6.0, 8.0 } });

    private static void AssertMatrix(GeneralMatrix actual, int rows, int columns, params double[] values)
    {
        Assert.AreEqual(rows, actual.RowDimension);
        Assert.AreEqual(columns, actual.ColumnDimension);
        Assert.AreEqual(rows * columns, values.Length);
        for (int row = 0; row < rows; row++)
        {
            for (int column = 0; column < columns; column++)
            {
                Assert.AreEqual(values[row * columns + column], actual.GetElement(row, column), 1e-12);
            }
        }
    }

    [TestMethod]
    public void PackedConstructorsRespectTheirStorageOrderAndOffset()
    {
        double[] values = { 1, 2, 3, 4, 5, 6 };
        var columnPacked = new GeneralMatrix(values, 2);
        var rowPacked = new GeneralMatrix(2, values);
        AssertMatrix(columnPacked, 2, 3, 1, 3, 5, 2, 4, 6);
        AssertMatrix(rowPacked, 2, 3, 1, 2, 3, 4, 5, 6);
        AssertMatrix(new GeneralMatrix(new[] { -99.0, 1, 2, 3, 4, 5, 6, 99 }, 2, 3, 1), 2, 3, 1, 3, 5, 2, 4, 6);
        values[0] = 100;
        AssertMatrix(columnPacked, 2, 3, 1, 3, 5, 2, 4, 6);
        AssertMatrix(rowPacked, 2, 3, 1, 2, 3, 4, 5, 6);
    }

    [TestMethod]
    public void PackedConstructorsRejectInsufficientOrNondivisibleInput()
    {
        Assert.ThrowsExactly<ArgumentException>(() => new GeneralMatrix(new[] { 1.0, 2, 3 }, 2));
        Assert.ThrowsExactly<ArgumentException>(() => new GeneralMatrix(2, new[] { 1.0, 2, 3 }));
        Assert.ThrowsExactly<ArgumentException>(() => new GeneralMatrix(new[] { 1.0, 2, 3 }, 2, 2, 0));
        Assert.ThrowsExactly<ArgumentException>(() => new GeneralMatrix(new[] { 1.0 }, 0));
        Assert.ThrowsExactly<ArgumentException>(() => new GeneralMatrix(0, new[] { 1.0 }));
        AssertMatrix(new GeneralMatrix(System.Array.Empty<double>(), 0), 0, 0);
        AssertMatrix(new GeneralMatrix(0, System.Array.Empty<double>()), 0, 0);
    }

    [TestMethod]
    public void ConstructorsSupportConstantAndEmptyMatrices()
    {
        AssertMatrix(new GeneralMatrix(2, 3), 2, 3, 0, 0, 0, 0, 0, 0);
        AssertMatrix(new GeneralMatrix(2, 3, -2), 2, 3, -2, -2, -2, -2, -2, -2);
        AssertMatrix(new GeneralMatrix(0, 3), 0, 3);
        AssertMatrix(new GeneralMatrix(3, 0), 3, 0);
        Assert.AreEqual(0.0, new GeneralMatrix(0, 3).Norm1());
        Assert.AreEqual(0.0, new GeneralMatrix(3, 0).NormInf());
        Assert.AreEqual(0.0, new GeneralMatrix(0, 0).NormF());
        Assert.AreEqual(0.0, new GeneralMatrix(0, 0).Trace());
    }

    [TestMethod]
    public void ArrayConstructorsAliasStorageWhileFactoriesAndCopiesAreIndependent()
    {
        double[][] storage = { new[] { 1.0, 2 }, new[] { 3.0, 4 } };
        var aliased = new GeneralMatrix(storage);
        var uncheckedShape = new GeneralMatrix(storage, 2, 2);
        var created = GeneralMatrix.Create(storage);
        var copied = aliased.Copy();
        var cloned = (GeneralMatrix)aliased.Clone();
        double[][] arrayCopy = aliased.ArrayCopy;
        double[] rowCopy = aliased.RowPackedCopy;
        double[] columnCopy = aliased.ColumnPackedCopy;
        Assert.AreSame(storage, aliased.Array);
        Assert.AreSame(storage, uncheckedShape.Array);
        storage[0][0] = 9;
        Assert.AreEqual(9.0, aliased.GetElement(0, 0));
        Assert.AreEqual(9.0, uncheckedShape.GetElement(0, 0));
        AssertMatrix(created, 2, 2, 1, 2, 3, 4);
        AssertMatrix(copied, 2, 2, 1, 2, 3, 4);
        AssertMatrix(cloned, 2, 2, 1, 2, 3, 4);
        Assert.AreEqual(1.0, arrayCopy[0][0]);
        CollectionAssert.AreEqual(new[] { 1.0, 2, 3, 4 }, rowCopy);
        CollectionAssert.AreEqual(new[] { 1.0, 3, 2, 4 }, columnCopy);
        arrayCopy[1][1] = -100;
        Assert.AreEqual(4.0, aliased.GetElement(1, 1));
    }

    [TestMethod]
    public void ArithmeticAndOperatorsPreserveInputsAndComputeExpectedValues()
    {
        var left = Sample();
        var right = new GeneralMatrix(2, 2, 2);
        AssertMatrix(left + right, 2, 2, 4, -2, 8, 10);
        AssertMatrix(left - right, 2, 2, 0, -6, 4, 6);
        AssertMatrix(left * right, 2, 2, -4, -4, 28, 28);
        AssertMatrix(left.UnaryMinus(), 2, 2, -2, 4, -6, -8);
        AssertMatrix(left.Multiply(-0.5), 2, 2, -1, 2, -3, -4);
        AssertMatrix(left.ArrayMultiply(right), 2, 2, 4, -8, 12, 16);
        AssertMatrix(left.ArrayRightDivide(right), 2, 2, 1, -2, 3, 4);
        AssertMatrix(left.ArrayLeftDivide(right), 2, 2, 1, -0.5, 1.0 / 3, 0.25);
        AssertMatrix(left, 2, 2, 2, -4, 6, 8);
        AssertMatrix(right, 2, 2, 2, 2, 2, 2);
    }

    [TestMethod]
    public void InPlaceArithmeticReturnsAndMutatesTheReceiver()
    {
        var source = Sample();
        var two = new GeneralMatrix(2, 2, 2);
        Assert.AreSame(source, source.AddEquals(two));
        AssertMatrix(source, 2, 2, 4, -2, 8, 10);
        Assert.AreSame(source, source.SubtractEquals(two));
        AssertMatrix(source, 2, 2, 2, -4, 6, 8);
        Assert.AreSame(source, source.ArrayMultiplyEquals(two));
        AssertMatrix(source, 2, 2, 4, -8, 12, 16);
        Assert.AreSame(source, source.ArrayRightDivideEquals(two));
        AssertMatrix(source, 2, 2, 2, -4, 6, 8);
        Assert.AreSame(source, source.ArrayLeftDivideEquals(two));
        AssertMatrix(source, 2, 2, 1, -0.5, 1.0 / 3, 0.25);
        Assert.AreSame(source, source.MultiplyEquals(3));
        AssertMatrix(source, 2, 2, 3, -1.5, 1, 0.75);
    }

    [TestMethod]
    [DataRow(1, 2)]
    [DataRow(2, 1)]
    public void ElementwiseOperationsRejectDifferentRowsOrColumnsBeforeMutating(int rows, int columns)
    {
        var source = Sample();
        var other = new GeneralMatrix(rows, columns);
        Func<GeneralMatrix>[] operations =
        {
            () => source.Add(other), () => source.AddEquals(other),
            () => source.Subtract(other), () => source.SubtractEquals(other),
            () => source.ArrayMultiply(other), () => source.ArrayMultiplyEquals(other),
            () => source.ArrayLeftDivide(other), () => source.ArrayLeftDivideEquals(other),
            () => source.ArrayRightDivide(other), () => source.ArrayRightDivideEquals(other)
        };
        foreach (var operation in operations)
        {
            Assert.ThrowsExactly<ArgumentException>(() => operation());
            AssertMatrix(source, 2, 2, 2, -4, 6, 8);
        }
        Assert.ThrowsExactly<ArgumentException>(() => source.Multiply(new GeneralMatrix(3, 1)));
    }

    [TestMethod]
    public void ElementwiseDivisionRetainsIeeeZeroAndInfinityBehavior()
    {
        var numerator = new GeneralMatrix(new[] { new[] { 0.0, 1, -1 } });
        var quotient = numerator.ArrayRightDivide(new GeneralMatrix(1, 3));
        Assert.IsTrue(double.IsNaN(quotient.GetElement(0, 0)));
        Assert.AreEqual(double.PositiveInfinity, quotient.GetElement(0, 1));
        Assert.AreEqual(double.NegativeInfinity, quotient.GetElement(0, 2));
    }

    [TestMethod]
    public void RectangularMultiplicationTransposeAndNormsAreCorrect()
    {
        var matrix = new GeneralMatrix(new[] { new[] { 1.0, -2, 3 }, new[] { -4.0, 5, -6 } });
        AssertMatrix(matrix.Transpose(), 3, 2, 1, -4, -2, 5, 3, -6);
        AssertMatrix(matrix.Multiply(matrix.Transpose()), 2, 2, 14, -32, -32, 77);
        Assert.AreEqual(9.0, matrix.Norm1());
        Assert.AreEqual(15.0, matrix.NormInf());
        Assert.AreEqual(Math.Sqrt(91), matrix.NormF(), 1e-12);
        Assert.AreEqual(6.0, matrix.Trace());
        AssertMatrix(GeneralMatrix.Identity(2, 3), 2, 3, 1, 0, 0, 0, 1, 0);
        AssertMatrix(GeneralMatrix.Identity(3, 2), 3, 2, 1, 0, 0, 1, 0, 0);
        AssertMatrix(new GeneralMatrix(2, 0).Multiply(new GeneralMatrix(0, 3)), 2, 3, 0, 0, 0, 0, 0, 0);
    }

    [TestMethod]
    public void IndexedSubmatricesHonorOrderingAndAreIndependent()
    {
        var matrix = new GeneralMatrix(3, new[] { 1.0, 2, 3, 4, 5, 6, 7, 8, 9 });
        AssertMatrix(matrix.GetMatrix(new[] { 2, 0 }, new[] { 2, 0 }), 2, 2, 9, 7, 3, 1);
        AssertMatrix(matrix.GetMatrix(1, 2, new[] { 2, 0 }), 2, 2, 6, 4, 9, 7);
        AssertMatrix(matrix.GetMatrix(new[] { 2, 0 }, 1, 2), 2, 2, 8, 9, 2, 3);
        var submatrix = matrix.GetMatrix(0, 1, 0, 1);
        AssertMatrix(submatrix, 2, 2, 1, 2, 4, 5);
        submatrix.SetElement(0, 0, -1);
        Assert.AreEqual(1.0, matrix.GetElement(0, 0));
        AssertMatrix(matrix.GetMatrix(System.Array.Empty<int>(), new[] { 0 }), 0, 1);
        AssertMatrix(matrix.GetMatrix(new[] { 0 }, System.Array.Empty<int>()), 1, 0);
    }

    [TestMethod]
    public void SetMatrixOverloadsWriteOnlySelectedElements()
    {
        var patch = new GeneralMatrix(2, new[] { 1.0, 2, 3, 4 });
        var matrix = new GeneralMatrix(3, 3);
        matrix.SetMatrix(1, 2, 1, 2, patch);
        AssertMatrix(matrix, 3, 3, 0, 0, 0, 0, 1, 2, 0, 3, 4);
        matrix = new GeneralMatrix(3, 3);
        matrix.SetMatrix(new[] { 2, 0 }, new[] { 2, 0 }, patch);
        AssertMatrix(matrix, 3, 3, 4, 0, 3, 0, 0, 0, 2, 0, 1);
        matrix = new GeneralMatrix(3, 3);
        matrix.SetMatrix(new[] { 2, 0 }, 1, 2, patch);
        AssertMatrix(matrix, 3, 3, 0, 3, 4, 0, 0, 0, 0, 1, 2);
        matrix = new GeneralMatrix(3, 3);
        matrix.SetMatrix(1, 2, new[] { 2, 0 }, patch);
        AssertMatrix(matrix, 3, 3, 0, 0, 0, 2, 0, 1, 4, 0, 3);
    }

    [TestMethod]
    public void SubmatrixFailuresRetainTheUnderlyingIndexException()
    {
        var matrix = Sample();
        Action[] actions =
        {
            () => matrix.GetMatrix(-1, 0, 0, 1),
            () => matrix.GetMatrix(new[] { 0 }, new[] { -1 }),
            () => matrix.GetMatrix(0, 1, new[] { -1 }),
            () => matrix.GetMatrix(new[] { -1 }, 0, 1),
            () => matrix.SetMatrix(-1, 0, 0, 1, Sample()),
            () => matrix.SetMatrix(new[] { -1 }, new[] { 0 }, Sample()),
            () => matrix.SetMatrix(new[] { -1 }, 0, 1, Sample()),
            () => matrix.SetMatrix(0, 1, new[] { -1 }, Sample())
        };
        foreach (var action in actions)
        {
            var exception = Assert.ThrowsExactly<IndexOutOfRangeException>(action);
            Assert.AreEqual("Submatrix indices", exception.Message);
            Assert.IsInstanceOfType<IndexOutOfRangeException>(exception.InnerException);
        }
        Assert.ThrowsExactly<IndexOutOfRangeException>(() => matrix.GetElement(-1, 0));
        Assert.ThrowsExactly<IndexOutOfRangeException>(() => matrix.SetElement(0, -1, 42));
    }

    [TestMethod]
    public void EqualityCoversValueReferenceTypeAndDimensionDifferences()
    {
        var source = Sample();
        var copy = source.Copy();
        Assert.IsTrue(source.Equals(copy));
        Assert.IsTrue(source.Equals((object)copy));
        Assert.IsTrue(source.Equals(source));
        Assert.IsTrue(source.Equals((object)source));
        Assert.IsTrue(source == copy);
        Assert.IsFalse(source != copy);
        Assert.IsFalse(source.Equals((GeneralMatrix)null!));
        Assert.IsFalse(source.Equals((object)null!));
        Assert.IsFalse(source.Equals("matrix"));
        Assert.IsFalse(source.Equals(new GeneralMatrix(1, 2)));
        Assert.IsFalse(source.Equals(new GeneralMatrix(2, 1)));
        copy.SetElement(0, 1, 4);
        Assert.IsFalse(source.Equals(copy));
        Assert.IsFalse(source == copy);
        Assert.IsTrue(source != copy);
        // Legacy hashes use storage identity and dimensions, so value mutations leave them unchanged.
        int originalHash = source.GetHashCode();
        source.SetElement(0, 0, -100);
        Assert.AreEqual(originalHash, source.GetHashCode());
    }

    [TestMethod]
    public void RandomOverloadsRespectRangesAndDegenerateBounds()
    {
        foreach (double value in GeneralMatrix.Random(4, 5).RowPackedCopy)
        {
            Assert.IsTrue(value >= 0 && value < 1);
        }
        foreach (double value in GeneralMatrix.Random(4, 5, -4.5, -1.5).RowPackedCopy)
        {
            Assert.IsTrue(value >= -4.5 && value < -1.5);
        }
        foreach (double value in GeneralMatrix.Random(4, 5, -4, 2).RowPackedCopy)
        {
            Assert.IsTrue(value >= -4 && value < 2);
            Assert.AreEqual(Math.Truncate(value), value);
        }
        AssertMatrix(GeneralMatrix.Random(1, 2, 3.5, 3.5), 1, 2, 3.5, 3.5);
        AssertMatrix(GeneralMatrix.Random(1, 2, 3, 3), 1, 2, 3, 3);
        Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GeneralMatrix.Random(1, 1, 3, 2));
        AssertMatrix(GeneralMatrix.Random(0, 2), 0, 2);
    }

    [TestMethod]
    public void SolveTransposePreservesLegacyTransposedResult()
    {
        var coefficients = new GeneralMatrix(new[] { new[] { 2.0, 1 }, new[] { 0.0, 3 } });
        var rightHandSide = new GeneralMatrix(new[] { new[] { 8.0, 19 } });
        var legacyResult = coefficients.SolveTranspose(rightHandSide);
        // The implementation returns X-transpose although its XML documentation says X.
        AssertMatrix(legacyResult, 2, 1, 4, 5);
        AssertMatrix(legacyResult.Transpose().Multiply(coefficients), 1, 2, 8, 19);
    }

    [TestMethod]
    public void DisposeIsIdempotentAndPreservesManagedMatrixValues()
    {
        var matrix = Sample();
        matrix.Dispose();
        matrix.Dispose();
        AssertMatrix(matrix, 2, 2, 2, -4, 6, 8);
        matrix.SetElement(0, 0, 3);
        Assert.AreEqual(3.0, matrix.GetElement(0, 0));
    }

    [TestMethod]
    public void LegacySerializationContractEmitsNoPayload()
    {
        // Characterize the existing empty ISerializable callback without invoking a formatter.
#pragma warning disable SYSLIB0050 // Preserve inspection of the legacy serialization contract.
        var info = new SerializationInfo(typeof(GeneralMatrix), new FormatterConverter());
        ((ISerializable)Sample()).GetObjectData(info, default);
        Assert.AreEqual(0, info.MemberCount);
        Assert.AreEqual(typeof(GeneralMatrix), info.ObjectType);
#pragma warning restore SYSLIB0050
    }
}
