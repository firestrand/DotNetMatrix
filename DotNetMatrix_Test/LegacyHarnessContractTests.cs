using System;
using System.Globalization;
using System.IO;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class LegacyHarnessContractTests
{
    [TestMethod]
    public void VectorCheckUsesScalarRelativeToleranceAndBothNearZeroCases()
    {
        double epsilon = Math.ScaleB(1, -52);
        TestMatrix.Check(new[] { 0.0, 5 * epsilon, 1 }, new[] { 5 * epsilon, 0.0, 1 + 5 * epsilon });
        TestMatrix.Check(Array.Empty<double>(), Array.Empty<double>());
        SystemException error = Assert.ThrowsExactly<SystemException>(() =>
            TestMatrix.Check(new[] { 1.0, 2 }, new[] { 1.0 + 20 * epsilon, 2 }));
        Assert.IsTrue(error.Message.Contains("difference x-y is too large", StringComparison.Ordinal));
        Assert.ThrowsExactly<SystemException>(() => TestMatrix.Check(new[] { 0.0 }, new[] { 20 * epsilon }));
    }

    [TestMethod]
    public void VectorCheckRejectsLengthMismatchBeforeComparingElements()
    {
        SystemException error = Assert.ThrowsExactly<SystemException>(() =>
            TestMatrix.Check(new[] { 1.0 }, new[] { 1.0, 2 }));
        Assert.AreEqual("Attempt to compare vectors of different lengths", error.Message);
    }

    [TestMethod]
    public void ArrayCheckUsesMatrixNormToleranceWithoutChangingInputs()
    {
        double epsilon = Math.ScaleB(1, -52);
        double[][] left = { new[] { 1.0, 2 }, new[] { 3.0, 4 } };
        double[][] close = { new[] { 1.0 + 100 * epsilon, 2 }, new[] { 3.0, 4 } };
        TestMatrix.Check(left, close);
        TestMatrix.Check(new[] { new[] { 0.0 } }, new[] { new[] { 5 * epsilon } });
        TestMatrix.Check(new[] { new[] { 5 * epsilon } }, new[] { new[] { 0.0 } });
        CollectionAssert.AreEqual(new[] { 1.0, 2, 3, 4 }, new GeneralMatrix(left).RowPackedCopy);
        CollectionAssert.AreEqual(new[] { 1.0 + 100 * epsilon, 2, 3, 4 }, new GeneralMatrix(close).RowPackedCopy);
        SystemException error = Assert.ThrowsExactly<SystemException>(() =>
            TestMatrix.Check(left, new[] { new[] { 1.01, 2 }, new[] { 3.0, 4 } }));
        Assert.IsTrue(error.Message.Contains("norm of (X-Y) is too large", StringComparison.Ordinal));
    }

    [TestMethod]
    public void ArrayCheckPreservesRaggedInputAndNonzeroShapeErrors()
    {
        Assert.ThrowsExactly<ArgumentException>(() =>
            TestMatrix.Check(new[] { new[] { 1.0 }, new[] { 2.0, 3 } }, new[] { new[] { 1.0 } }));
        Assert.ThrowsExactly<ArgumentException>(() =>
            TestMatrix.Check(new[] { new[] { 1.0, 2 } }, new[] { new[] { 1.0 }, new[] { 2.0 } }));
    }

    [TestMethod]
    public void WarningReportingIncrementsItsCounterAndRetainsTheLegacyMessage()
    {
        TextWriter original = Console.Out;
        using var output = new StringWriter(CultureInfo.InvariantCulture);
        try
        {
            Console.SetOut(output);
            int count = TestMatrix.TryWarning(4, "factorization", "tolerance warning");
            Assert.AreEqual(5, count);
            Assert.AreEqual(">    factorization*** warning ***\n>      Message: tolerance warning\n", output.ToString());
        }
        finally
        {
            Console.SetOut(original);
        }
    }
}
