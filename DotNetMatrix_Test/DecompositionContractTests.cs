using System;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class DecompositionContractTests
{
    [TestMethod]
    public void DecompositionsRejectNullEmptyAndNonfiniteInputs()
    {
        Action<GeneralMatrix>[] constructors =
        [
            a => _ = new CholeskyDecomposition(a), a => _ = new LUDecomposition(a),
            a => _ = new QRDecomposition(a), a => _ = new SingularValueDecomposition(a),
            a => _ = new EigenvalueDecomposition(a)
        ];
        foreach (var construct in constructors)
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => construct(null!));
            foreach (var empty in new[] { new GeneralMatrix(0, 0), new GeneralMatrix(0, 2), new GeneralMatrix(2, 0) })
                Assert.ThrowsExactly<ArgumentException>(() => construct(empty));
            foreach (double invalid in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
                Assert.ThrowsExactly<ArgumentException>(() => construct(new GeneralMatrix(new[] { new[] { invalid } })));
            Assert.ThrowsExactly<ArgumentException>(() => construct(new GeneralMatrix(new[] { new[] { 1d, 2d } }, 1, 1)));
            Assert.ThrowsExactly<ArgumentException>(() => construct(new GeneralMatrix(new double[0][], 1, 1)));
            Assert.ThrowsExactly<ArgumentException>(() => construct(new GeneralMatrix(new[] { (double[])null! }, 1, 1)));
        }
    }

    [TestMethod]
    public void RectangularDomainsAndTallLuSolveAreExplicit()
    {
        var wide = Matrix([1, 0, 2], [0, 1, 3]);
        var tall = wide.Transpose();
        Assert.ThrowsExactly<ArgumentException>(() => new LUDecomposition(wide));
        Assert.ThrowsExactly<ArgumentException>(() => new QRDecomposition(wide));
        Assert.ThrowsExactly<ArgumentException>(() => new CholeskyDecomposition(tall));
        Assert.ThrowsExactly<ArgumentException>(() => new EigenvalueDecomposition(tall));
        Assert.ThrowsExactly<ArgumentException>(() => new EigenvalueDecomposition(wide));
        Assert.ThrowsExactly<ArgumentException>(() => new LUDecomposition(tall).Solve(new GeneralMatrix(3, 1)));
        Assert.AreEqual(3, new QRDecomposition(tall).Q.RowDimension);
    }

    [TestMethod]
    public void SolversValidateRightHandSideAndSupportNoColumns()
    {
        var identity = GeneralMatrix.Identity(2, 2);
        Func<GeneralMatrix, GeneralMatrix>[] solves =
        [new CholeskyDecomposition(identity).Solve, new LUDecomposition(identity).Solve, new QRDecomposition(identity).Solve];
        foreach (var solve in solves)
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => solve(null!));
            Assert.ThrowsExactly<ArgumentException>(() => solve(Matrix([double.NaN], [0])));
            var emptyRight = solve(new GeneralMatrix(2, 0));
            Assert.AreEqual(2, emptyRight.RowDimension);
            Assert.AreEqual(0, emptyRight.ColumnDimension);
        }
    }

    [TestMethod]
    public void IterationBudgetsAreValidatedAndBothEigenAlgorithmsAreBounded()
    {
        var symmetric = Matrix([2, 1, 0], [1, 2, 1], [0, 1, 2]);
        var nonsymmetric = Matrix([1, 2, 3], [4, 5, 6], [7, 8, 10]);
        foreach (int invalid in new[] { 0, -1 })
        {
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new SingularValueDecomposition(symmetric, invalid));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new EigenvalueDecomposition(symmetric, invalid));
        }
        foreach (var input in new[] { symmetric, nonsymmetric })
        {
            var error = Assert.ThrowsExactly<DecompositionConvergenceException>(() => new EigenvalueDecomposition(input, 1));
            Assert.AreEqual(nameof(EigenvalueDecomposition), error.Algorithm);
            Assert.AreEqual(1, error.MaxIterations);
            var eig = new EigenvalueDecomposition(input, 1000);
            Near(input.Multiply(eig.GetV()), eig.GetV().Multiply(eig.D));
        }
        var svdError = Assert.ThrowsExactly<DecompositionConvergenceException>(() => new SingularValueDecomposition(nonsymmetric, 1));
        Assert.AreEqual(nameof(SingularValueDecomposition), svdError.Algorithm);
        Assert.AreEqual(1, svdError.MaxIterations);
        Assert.ThrowsExactly<DecompositionConvergenceException>(() => new SingularValueDecomposition(nonsymmetric.GetMatrix(0, 1, 0, 2), 1));
    }

    [TestMethod]
    public void WideDenseAndIllConditionedSvdReconstructsAndKeepsInputs()
    {
        for (int rows = 1; rows <= 6; rows++)
        {
            var input = new GeneralMatrix(rows, rows + 3);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < input.ColumnDimension; j++)
                    input.SetElement(i, j, Math.Sin((i + 1) * (j + 2)) + (i == j ? 1e-6 : 0));
            var copy = input.Copy();
            var svd = new SingularValueDecomposition(input, 1000);
            Assert.AreEqual(rows, svd.GetU().ColumnDimension);
            Assert.AreEqual(rows, svd.GetV().ColumnDimension);
            Assert.AreEqual(rows, svd.SingularValues.Length);
            Near(input, svd.GetU().Multiply(svd.S).Multiply(svd.GetV().Transpose()));
            Near(GeneralMatrix.Identity(rows, rows), svd.GetU().Transpose().Multiply(svd.GetU()));
            Near(GeneralMatrix.Identity(rows, rows), svd.GetV().Transpose().Multiply(svd.GetV()));
            Near(copy, input);
        }
        var illConditioned = Matrix([1, 1, 1], [1, 1 + 1e-10, 1], [1, 1, 1 + 1e-12]);
        var factors = new SingularValueDecomposition(illConditioned, 1000);
        Near(illConditioned, factors.GetU().Multiply(factors.S).Multiply(factors.GetV().Transpose()));
    }

    [TestMethod]
    public void CholeskyAndSvdCopyAccessorsDoNotChangeLegacyAliases()
    {
        var chol = new CholeskyDecomposition(Matrix([4, 2], [2, 10]));
        chol.GetLCopy().SetElement(0, 0, -100);
        Assert.AreEqual(2d, chol.GetL().GetElement(0, 0));
        chol.GetL().SetElement(0, 0, 3);
        Assert.AreEqual(3d, chol.GetL().GetElement(0, 0));
        foreach (var input in new[] { Matrix([3, 0], [0, 2]), Matrix([3, 0, 0], [0, 2, 0]) })
        {
            var svd = new SingularValueDecomposition(input);
            svd.GetUCopy().SetElement(0, 0, -100);
            svd.GetVCopy().SetElement(0, 0, -100);
            svd.GetSingularValuesCopy()[0] = -100;
            Assert.AreNotEqual(-100d, svd.GetU().GetElement(0, 0));
            Assert.AreNotEqual(-100d, svd.GetV().GetElement(0, 0));
            Assert.AreEqual(3d, svd.SingularValues[0]);
            svd.GetU().SetElement(0, 0, -9);
            svd.GetV().SetElement(0, 0, -8);
            svd.SingularValues[0] = 7;
            Assert.AreEqual(-9d, svd.GetU().GetElement(0, 0));
            Assert.AreEqual(-8d, svd.GetV().GetElement(0, 0));
            Assert.AreEqual(7d, svd.SingularValues[0]);
            svd.S.SetElement(0, 0, -10);
            Assert.AreEqual(7d, svd.SingularValues[0]);
        }
    }

    [TestMethod]
    public void EigenCopyAccessorsPreserveIndependentCopiesAndLegacyAliases()
    {
        var eig = new EigenvalueDecomposition(Matrix([0, -1], [1, 0]));
        eig.GetVCopy().SetElement(0, 0, -100);
        eig.GetRealEigenvaluesCopy()[0] = -100;
        eig.GetImagEigenvaluesCopy()[0] = -100;
        Assert.AreNotEqual(-100d, eig.GetV().GetElement(0, 0));
        Assert.AreEqual(0d, eig.RealEigenvalues[0]);
        Assert.AreEqual(1d, eig.ImagEigenvalues[0]);
        eig.GetV().SetElement(0, 0, -9);
        eig.RealEigenvalues[0] = 7;
        eig.ImagEigenvalues[0] = 8;
        Assert.AreEqual(-9d, eig.GetV().GetElement(0, 0));
        Assert.AreEqual(7d, eig.RealEigenvalues[0]);
        Assert.AreEqual(8d, eig.ImagEigenvalues[0]);
        eig.D.SetElement(0, 0, -100);
        Assert.AreEqual(7d, eig.RealEigenvalues[0]);
    }

    [TestMethod]
    public void LuAndQrFactorAccessorsProduceIndependentMatrices()
    {
        var input = Matrix([2, 1], [1, 2], [0, 1]);
        var lu = new LUDecomposition(input);
        lu.L.SetElement(0, 0, -100);
        lu.U.SetElement(0, 0, -100);
        lu.DoublePivot[0] = -100;
        Assert.AreNotEqual(-100d, lu.L.GetElement(0, 0));
        Assert.AreNotEqual(-100d, lu.U.GetElement(0, 0));
        Assert.AreNotEqual(-100d, lu.DoublePivot[0]);
        var qr = new QRDecomposition(input);
        qr.Q.SetElement(0, 0, -100);
        qr.R.SetElement(0, 0, -100);
        qr.H.SetElement(0, 0, -100);
        Assert.AreNotEqual(-100d, qr.Q.GetElement(0, 0));
        Assert.AreNotEqual(-100d, qr.R.GetElement(0, 0));
        Assert.AreNotEqual(-100d, qr.H.GetElement(0, 0));
        Near(input, qr.Q.Multiply(qr.R));
    }

    private static GeneralMatrix Matrix(params double[][] rows) => new(rows);

    private static void Near(GeneralMatrix expected, GeneralMatrix actual)
    {
        Assert.AreEqual(expected.RowDimension, actual.RowDimension);
        Assert.AreEqual(expected.ColumnDimension, actual.ColumnDimension);
        Assert.IsTrue(expected.Subtract(actual).NormF() <= 1e-9 * Math.Max(1, expected.NormF()));
    }
}
