using System.Collections.Generic;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class MatrixContractTests
{
    private sealed class DerivedMatrix() : GeneralMatrix(1, 1);

    [TestMethod]
    public void EqualityIsSymmetricForNullTypesNaNAndSignedZero()
    {
        GeneralMatrix? absent = null;
        var a = new GeneralMatrix(new[] { new[] { double.NaN, -0.0 } });
        var b = new GeneralMatrix(new[] { new[] { double.NaN, 0.0 } });
        Assert.IsTrue(absent! == null!);
        Assert.IsFalse(absent! != null!);
        Assert.IsFalse(absent! == a);
        Assert.IsFalse(a == absent!);
        Assert.IsTrue(absent! != a);
        Assert.IsTrue(a != absent!);
        Assert.IsTrue(a.Equals(b));
        Assert.IsTrue(a.Equals((object)b));
        Assert.AreEqual(a.GetHashCode(), b.GetHashCode());
        var dictionary = new Dictionary<GeneralMatrix, string> { [a] = "value" };
        Assert.AreEqual("value", dictionary[b]);
        var plain = new GeneralMatrix(1, 1);
        var derived = new DerivedMatrix();
        Assert.IsFalse(plain.Equals(derived));
        Assert.IsFalse(derived.Equals(plain));
        Assert.IsTrue(derived.Equals(new DerivedMatrix()));
        Assert.IsTrue(derived.Equals((object)new DerivedMatrix()));
    }

    [TestMethod]
    public void RectangularTransposeSolveSatisfiesLeastSquaresOptimality()
    {
        var a = new GeneralMatrix(new[] { new[] { 1.0, 0, 0 }, new[] { 0.0, 1, 0 } });
        var b = new GeneralMatrix(new[] { new[] { 2.0, 3, 4 } });
        var x = a.SolveTranspose(b);
        Assert.AreEqual(1, x.RowDimension);
        Assert.AreEqual(2, x.ColumnDimension);
        Assert.AreEqual(2, x.GetElement(0, 0), 1e-12);
        Assert.AreEqual(3, x.GetElement(0, 1), 1e-12);
        var residual = x.Multiply(a).Subtract(b);
        Assert.AreEqual(4, residual.NormF(), 1e-12);
        Assert.AreEqual(0, residual.Multiply(a.Transpose()).NormF(), 1e-12);
        Assert.AreEqual(4, b.GetElement(0, 2));
        Assert.AreEqual(1, a.GetElement(0, 0));
    }
}
