using System;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class CopyContractTests
{
    [TestMethod]
    public void CopyPreservesFloatingPointBitsAndOwnsEveryRow()
    {
        double[] row = { -0.0, double.PositiveInfinity, double.NegativeInfinity,
            BitConverter.Int64BitsToDouble(0x7ff8000000000042),
            BitConverter.Int64BitsToDouble(0x7ff0000000000001), double.Epsilon, double.MaxValue };
        Array.Resize(ref row, 16);
        var source = new GeneralMatrix(new[] { row, row });
        GeneralMatrix copy = source.Copy();
        Assert.AreEqual(2, copy.RowDimension);
        Assert.AreEqual(row.Length, copy.ColumnDimension);
        Assert.AreNotSame(source.Array, copy.Array);
        Assert.AreNotSame(copy.Array[0], copy.Array[1]);
        for (int i = 0; i < 2; i++)
        {
            Assert.AreNotSame(row, copy.Array[i]);
            for (int j = 0; j < row.Length; j++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(row[j]), BitConverter.DoubleToInt64Bits(copy.Array[i][j]));
        }
        copy.Array[0][0] = 123;
        Assert.AreEqual(long.MinValue, BitConverter.DoubleToInt64Bits(row[0]));
        Assert.AreEqual(long.MinValue, BitConverter.DoubleToInt64Bits(copy.Array[1][0]));
    }

    [TestMethod]
    public void CopyIgnoresStorageOutsideLogicalDimensions()
    {
        var source = new GeneralMatrix(new[] { new[] { 1.0, 2, 99 }, new[] { 3.0, 4, 99 }, null! }, 2, 2);
        CollectionAssert.AreEqual(new[] { 1.0, 2, 3, 4 }, source.Copy().RowPackedCopy);
    }

    [TestMethod]
    public void BulkCopyIgnoresStorageBeyondTheLogicalRow()
    {
        double[] row = new double[17];
        for (int i = 0; i < 16; i++) row[i] = i + 0.5;
        row[16] = double.NaN;
        var source = new GeneralMatrix(new[] { row }, 1, 16);
        GeneralMatrix copy = source.Copy();
        Assert.AreEqual(16, copy.Array[0].Length);
        for (int i = 0; i < 16; i++) Assert.AreEqual(i + 0.5, copy.Array[0][i]);
    }

    [TestMethod]
    [DataRow(2)]
    [DataRow(16)]
    public void CopyRetainsShortRowException(int columns)
    {
        var source = new GeneralMatrix(new[] { new[] { 1.0 } }, 1, columns);
        Assert.ThrowsExactly<IndexOutOfRangeException>(() => source.Copy());
    }

    [TestMethod]
    [DataRow(2)]
    [DataRow(16)]
    public void CopyRetainsMissingRowException(int columns)
    {
        var source = new GeneralMatrix(Array.Empty<double[]>(), 1, columns);
        Assert.ThrowsExactly<IndexOutOfRangeException>(() => source.Copy());
    }

    [TestMethod]
    [DataRow(2)]
    [DataRow(16)]
    public void CopyRetainsNullRowException(int columns)
    {
        var source = new GeneralMatrix(new double[][] { null! }, 1, columns);
        Assert.ThrowsExactly<NullReferenceException>(() => source.Copy());
    }

    [TestMethod]
    public void EmptyCopyDoesNotInspectIgnoredBorrowedStorage()
    {
        GeneralMatrix copy = new GeneralMatrix(null!, 3, 0).Copy();
        Assert.AreEqual(3, copy.RowDimension);
        Assert.AreEqual(0, copy.ColumnDimension);
        Assert.AreEqual(3, copy.Array.Length);
        foreach (double[] row in copy.Array) Assert.AreEqual(0, row.Length);
        GeneralMatrix noRows = new GeneralMatrix(null!, 0, 16).Copy();
        Assert.AreEqual(0, noRows.RowDimension);
        Assert.AreEqual(16, noRows.ColumnDimension);
    }

    [TestMethod]
    public void SingleColumnCopyRetainsBitsOwnershipAndExceptions()
    {
        var source = new GeneralMatrix(new[] { new[] { -0.0 }, new[] { double.Epsilon } });
        GeneralMatrix copy = source.Copy();
        Assert.AreEqual(long.MinValue, BitConverter.DoubleToInt64Bits(copy.Array[0][0]));
        Assert.AreEqual(double.Epsilon, copy.Array[1][0]);
        Assert.AreNotSame(source.Array[0], copy.Array[0]);
        Assert.ThrowsExactly<IndexOutOfRangeException>(() => new GeneralMatrix(new[] { Array.Empty<double>() }, 1, 1).Copy());
        Assert.ThrowsExactly<NullReferenceException>(() => new GeneralMatrix(new double[][] { null! }, 1, 1).Copy());
    }
}
