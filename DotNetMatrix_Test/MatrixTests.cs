using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test
{
    [TestClass]
    public class MatrixTests
    {
        GeneralMatrix A = null!, B = null!, C = null!, Z = null!, O = null!, R = null!, S = null!, M = null!, twos = null!;
        double[] columnwise { get { return new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 }; } }
        double[] rowwise { get { return new[] { 1.0, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 3.0, 6.0, 9.0, 12.0 }; } }
        double[][] avals { get { return new[] { new[] { 1.0, 4.0, 7.0, 10.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } }; } }

        //double[][] rankdef = avals;
        double[][] tvals { get { return new[] { new[] { 1.0, 2.0, 3.0 }, new[] { 4.0, 5.0, 6.0 }, new[] { 7.0, 8.0, 9.0 }, new[] { 10.0, 11.0, 12.0 } }; } }
        double[][] subavals { get { return new[] { new[] { 5.0, 8.0, 11.0 }, new[] { 6.0, 9.0, 12.0 } }; } }
        double[][] rvals { get { return new[] { new[] { 1.0, 4.0, 7.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } }; } }
        double[][] pvals { get { return new[] { new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 2.0, 3.0 }, new[] { 1.0, 3.0, 6.0 } }; } }
        double[][] ivals { get { return new[] { new[] { 1.0, 0.0, 0.0, 0.0 }, new[] { 0.0, 1.0, 0.0, 0.0 }, new[] { 0.0, 0.0, 1.0, 0.0 } }; } }
        double[][] evals = { new[] { 0.0, 1.0, 0.0, 0.0 }, new[] { 1.0, 0.0, 2e-7, 0.0 }, new[] { 0.0, -2e-7, 0.0, 1.0 }, new[] { 0.0, 0.0, 1.0, 0.0 } };
        double[][] square = { new[] { 166.0, 188.0, 210.0 }, new[] { 188.0, 214.0, 240.0 }, new[] { 210.0, 240.0, 270.0 } };
        double[][] sqSolution = { new[] { 13.0 }, new[] { 15.0 } };
        double[][] condmat = { new[] { 1.0, 3.0 }, new[] { 7.0, 9.0 } };
        int rows = 3, cols = 4;
        int nonconformld = 4; /* leading dimension which is valid, but nonconforming */
        int ib = 1, ie = 2, jb = 1, je = 3; /* index ranges for sub GeneralMatrix */
        int[] rowindexset = new int[] { 1, 2 };
        int[] badrowindexset = new int[] { 1, 3 };
        int[] columnindexset = new int[] { 1, 2, 3 };
        int[] badcolumnindexset = new int[] { 1, 2, 4 };

        //General
        /// <summary>Check norm of difference of Matrices. *</summary>

        private static bool check(GeneralMatrix x, GeneralMatrix y)
        {
            double eps = Math.Pow(2.0, -52.0);
            if (x.Norm1() == 0.0 & y.Norm1() < 10 * eps)
                return true;
            if (y.Norm1() == 0.0 & x.Norm1() < 10 * eps)
                return true;
            return x.Subtract(y).Norm1() <= 1000 * eps * Math.Max(x.Norm1(), y.Norm1());
        }
        [TestInitialize]
        public void InitializeArrays()
        {
            B = new GeneralMatrix(avals);
            M = new GeneralMatrix(2, 3, 0.0);
            S = new GeneralMatrix(columnwise, nonconformld);
            R = GeneralMatrix.Random(rows, cols);
            Z = new GeneralMatrix(rows, cols);
            A = R.Copy();
            C = A.Subtract(B);
            O = new GeneralMatrix(rows, cols, 1.0);
            twos = new GeneralMatrix(rows, cols, 2.0);
        }
        [TestMethod]
        public void RandomWithHighAndLowParametersGeneratesARandomMatrixWithValuesBetweenTheHighAndLowParameters()
        {
            double low = 2.0;
            double high = 5.0;
            var matrix = GeneralMatrix.Random(rows, cols, 2.0, 5.0);
            foreach (var val in matrix.ColumnPackedCopy)
            {
                Assert.IsTrue(val >= low && val <= high);
            }
        }
        /// <summary>
        /// check that exception is thrown in packed constructor with invalid length
        /// </summary>
        [TestMethod]
        public void ArgumentExceptionIsThrownInPackedConstructorWithInvalidLength()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                int invalid = 5; /* should trigger bad shape for construction with val */
                var a = new GeneralMatrix(columnwise, invalid);
                Assert.Fail("Expected construction to throw, but it returned a matrix.");
            });
        }
        /// <summary>
        /// check that exception is thrown in default constructor if input array is 'ragged' *
        /// </summary>
        [TestMethod]
        public void ArgumentExceptionIsThrownInConstructorWithRaggedInputArray()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                double[][] array = { new[] { 1.0, 4.0, 7.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } };
                var a = new GeneralMatrix(array);
                Assert.Fail("Expected construction to throw, but it returned a matrix.");
            });
        }
        /// <summary>
        /// check that exception is thrown in Create if input array is 'ragged' *
        /// </summary>
        [TestMethod]
        public void ArgumentExceptionIsThrownInCreateWithRaggedInputArray()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                double[][] array = { new[] { 1.0, 4.0, 7.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } };
                var a = GeneralMatrix.Create(array);
                Assert.Fail("Expected construction to throw, but it returned a matrix.");
            });
        }

        [TestMethod]
        public void CanBeConstructedFromArrayAndRowOffset()
        {
            double[] packedColumns = new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            int numRows = 3; /* leading dimension of intended test Matrices */
            var a = new GeneralMatrix(packedColumns, numRows);
            Assert.IsTrue(a.RowDimension == numRows);
            Assert.IsTrue(a.ColumnDimension == packedColumns.Length / numRows);
        }

        [TestMethod]
        public void ConstructionFromAPackedArrayDoesNotAltersOriginalArrayValues()
        {
            double[] packedColumns = new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            int numRows = 3; /* leading dimension of intended test Matrices */
            var b = new GeneralMatrix(packedColumns, numRows); //B still references array underneath
            var temp = b.GetElement(0, 0);
            packedColumns[0] = 0.0;
            Assert.IsTrue((temp - b.GetElement(0, 0)) == 0.0);
        }

        [TestMethod]
        public void ConstructionFromAnArrayAltersOriginalArrayValues()
        {
            double[][] array = avals;
            var b = new GeneralMatrix(array); //B still references avals underneath
            var temp = b.GetElement(0, 0);
            array[0][0] = 0.0;
            Assert.IsTrue((temp - b.GetElement(0, 0)) != 0.0);
        }

        [TestMethod]
        public void CreateFromAnArrayDoesNotAlterOriginalArrayValues()
        {
            double[][] array = avals;
            var b = GeneralMatrix.Create(array);
            var temp = b.GetElement(0, 0);
            array[0][0] = 0.0;
            Assert.IsTrue((temp - b.GetElement(0, 0)) == 0.0);
        }

        [TestMethod]
        public void CanCreateAnIdentityMatrixThroughIdentityMethod()
        {
            double[][] idVals = { new[] { 1.0, 0.0, 0.0, 0.0 }, new[] { 0.0, 1.0, 0.0, 0.0 }, new[] { 0.0, 0.0, 1.0, 0.0 } };
            var expected = new GeneralMatrix(idVals);
            var actual = GeneralMatrix.Identity(3, 4);
            Assert.IsTrue(expected == actual);
        }
        [TestMethod]
        public void CanAccessRowDimensions()
        {
            var matrix = new GeneralMatrix(avals);
            Assert.AreEqual(3, matrix.RowDimension);
        }
        [TestMethod]
        public void CanAccessColumnDimensions()
        {
            var matrix = new GeneralMatrix(avals);
            Assert.AreEqual(4, matrix.ColumnDimension);
        }
        [TestMethod]
        public void ArrayPropertyReturnsReferenceToUnderlyingArray()
        {
            var expected = avals;
            var matrix = new GeneralMatrix(expected);
            Assert.AreSame(expected, matrix.Array);
        }
        [TestMethod]
        public void ArrayCopyPropertyReturnsADeepCopyOfUnderlyingArray()
        {
            var matrix = new GeneralMatrix(avals);
            var copy = matrix.ArrayCopy;
            Assert.AreNotSame(matrix.Array, copy);
            for (int row = 0; row < matrix.RowDimension; row++)
            {
                Assert.AreNotSame(matrix.Array[row], copy[row]);
                Assert.IsTrue(matrix.Array[row].SequenceEqual(copy[row]));
            }
        }
        [TestMethod]
        public void ColumnPackedCopyPropertyReturnsAColumnPackedCopyOfUnderlyingArray()
        {
            var expected = new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            var matrix = new GeneralMatrix(avals);
            Assert.IsTrue(expected.SequenceEqual(matrix.ColumnPackedCopy));
            matrix = new GeneralMatrix(expected, 3);
            Assert.IsTrue(expected.SequenceEqual(matrix.ColumnPackedCopy));
            Assert.AreNotSame<object>(expected, matrix);
        }
        [TestMethod]
        public void RowPackedCopyPropertyReturnsARowPackedCopyOfUnderlyingArray()
        {
            var expected = new[] { 1.0, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 3.0, 6.0, 9.0, 12.0 };
            var matrix = new GeneralMatrix(avals);
            Assert.IsTrue(expected.SequenceEqual(matrix.RowPackedCopy));
            matrix = new GeneralMatrix(3, expected);
            Assert.IsTrue(expected.SequenceEqual(matrix.RowPackedCopy));
            Assert.AreNotSame<object>(expected, matrix);
        }
        [TestMethod]
        public void ZeroElementAccessThrowsAnIndexOutOfRangeExceptionWhenColumnValueIsEqualToColumnDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetElement(matrix.RowDimension - 1, matrix.ColumnDimension);
            });
        }
        [TestMethod]
        public void ZeroElementAccessThrowsAnIndexOutOfRangeExceptionWhenRowValueIsEqualToRowDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetElement(matrix.RowDimension, matrix.ColumnDimension - 1);
            });
        }
        [TestMethod]
        public void ZeroBasedGetElementMethodReturnsCorrectElement()
        {
            var matrix = new GeneralMatrix(avals);
            var actual = matrix.GetElement(matrix.RowDimension - 1, matrix.ColumnDimension - 1);
            Assert.AreEqual(12.0, actual);
        }
        [TestMethod]
        public void GetSubMatrixThrowsAnIndexOutOfRangeExceptionWhenFinalRowIndexExceedsRowDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(ib, ie + matrix.RowDimension + 1, jb, je);
            });
        }
        [TestMethod]
        public void GetSubMatrixThrowsAnIndexOutOfRangeExceptionWhenFinalColumnIndexExceedsColumnDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(ib, ie, jb, je + matrix.ColumnDimension + 1);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingColumnAndRowStartStopGetsCorrectSubMatrix()
        {
            var matrix = new GeneralMatrix(avals);
            var expected = new GeneralMatrix(subavals);
            var actual = matrix.GetMatrix(ib, ie, jb, je);
            Assert.AreEqual(expected, actual);
        }
        [TestMethod]
        public void GetSubMatrixThrowsAnIndexOutOfRangeExceptionWhenAColumnIndexInTheIndexArrayExceedsColumnDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(ib, ie, badcolumnindexset);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingColumnIndexArrayThrowsAnIndexOutOfRangeExceptionWhenAFinalRowIndexExceedsRowDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(ib, ie + matrix.RowDimension + 1, columnindexset);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingRowStartStopAndColumnIdexArrayGetsCorrectSubMatrix()
        {
            var matrix = new GeneralMatrix(avals);
            var expected = new GeneralMatrix(subavals);
            var actual = matrix.GetMatrix(ib, ie, columnindexset);
            Assert.AreEqual(expected, actual);
        }
        [TestMethod]
        public void GetSubMatrixThrowsAnIndexOutOfRangeExceptionWhenARowIndexInTheIndexArrayExceedsRowDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(badrowindexset, jb, je);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingRowIndexArrayThrowsAnIndexOutOfRangeExceptionWhenAFinalColumnIndexExceedsColumnDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(rowindexset, jb, je + matrix.ColumnDimension + 1);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingColumnStartStopAndRowIdexArrayGetsCorrectSubMatrix()
        {
            var matrix = new GeneralMatrix(avals);
            var expected = new GeneralMatrix(subavals);
            var actual = matrix.GetMatrix(rowindexset, jb, je);
            Assert.AreEqual(expected, actual);
        }
        [TestMethod]
        public void GetSubMatrixUsingRowAndColumnArraysThrowsAnIndexOutOfRangeExceptionWhenARowIndexInTheIndexArrayExceedsRowDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(badrowindexset, columnindexset);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingRowAndColumnArraysThrowsAnIndexOutOfRangeExceptionWhenAColumnIndexInTheIndexArrayExceedsColumnDimension()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.GetMatrix(rowindexset, badcolumnindexset);
            });
        }
        [TestMethod]
        public void GetSubMatrixUsingColumnAndRowIdexArrayGetsCorrectSubMatrix()
        {
            var matrix = new GeneralMatrix(avals);
            var expected = new GeneralMatrix(subavals);
            var actual = matrix.GetMatrix(rowindexset, columnindexset);
            Assert.AreEqual(expected, actual);
        }
        [TestMethod]
        public void SetElementWithInvalidRowIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.SetElement(matrix.RowDimension, matrix.ColumnDimension + 1, 0.0);
            });
        }
        [TestMethod]
        public void SetElementWithInvalidColumnIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                var matrix = new GeneralMatrix(avals);
                matrix.SetElement(matrix.RowDimension + 1, matrix.ColumnDimension, 0.0);
            });
        }
        [TestMethod]
        public void SetElementCorrectlySetsAMatrixElement()
        {
            var matrix = new GeneralMatrix(avals);
            matrix.SetElement(ib, jb, 0.0);
            Assert.AreEqual(0.0, matrix.GetElement(ib, jb));
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndInvalidRowStopIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(ib, ie + B.RowDimension + 1, jb, je, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndInvalidColumnStopIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(ib, ie, jb, je + B.ColumnDimension + 1, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixCorrectlySetsTheMatrixValues()
        {
            B.SetMatrix(ib, ie, jb, je, M);
            Assert.IsTrue(check(M.Subtract(B.GetMatrix(ib, ie, jb, je)), M));
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndColumnArrayAndInvalidRowStopIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(ib, ie + B.RowDimension + 1, columnindexset, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndInvalidColumnArrayIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(ib, ie, badcolumnindexset, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndColumnSetArrayCorrectlySetsTheMatrixValues()
        {
            B.SetMatrix(ib, ie, columnindexset, M);
            Assert.IsTrue(check(M.Subtract(B.GetMatrix(ib, ie, columnindexset)), M));
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndRowArrayAndInvalidColumnStopIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(rowindexset, jb, je + B.ColumnDimension + 1, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndInvalidRowArrayIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(badrowindexset, jb, je, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndRowSetArrayCorrectlySetsTheMatrixValues()
        {
            B.SetMatrix(rowindexset, jb, je, M);
            Assert.IsTrue(check(M.Subtract(B.GetMatrix(rowindexset, jb, je)), M));
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndRowArrayAndInvalidColumnArrayIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(rowindexset, badcolumnindexset, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixAndColumnArrayAndInvalidRowArrayIndexThrowsAnException()
        {
            Assert.ThrowsExactly<IndexOutOfRangeException>(() =>
            {
                B.SetMatrix(badrowindexset, columnindexset, M);
            });
        }
        [TestMethod]
        public void SetMatrixWithPassedInMatrixRowAndColumnSetArrayCorrectlySetsTheMatrixValues()
        {
            B.SetMatrix(rowindexset, columnindexset, M);
            Assert.IsTrue(check(M.Subtract(B.GetMatrix(rowindexset, columnindexset)), M));
        }
        [TestMethod]
        public void SubtractWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                S = A.Subtract(S);
            });
        }
        [TestMethod]
        public void SubtractMatrixFromItselfShouldResultInZeroMatrix()
        {
            var actual = A.Subtract(R).Norm1();
            Assert.AreEqual(0.0, actual);
        }
        [TestMethod]
        public void SubtractEqualsWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                A.SubtractEquals(S);
            });
        }
        [TestMethod]
        public void SubtractEqualsMatrixFromItselfShouldResultInZeroMatrix()
        {
            A.SubtractEquals(R);
            var actual = A.Norm1();
            Assert.AreEqual(0.0, actual);
        }
        [TestMethod]
        public void AddWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                S = A.Add(S);
            });
        }
        [TestMethod]
        public void AddMatrixGivesCorrectResult()
        {
            var actual = O.Copy().Add(O);
            Assert.AreEqual(twos, actual);
        }
        [TestMethod]
        public void AddEqualsWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                A.AddEquals(S);
            });
        }
        [TestMethod]
        public void AddEqualsMatrixGivesCorrectResult()
        {
            var actual = O.Copy();
            actual.AddEquals(O);
            Assert.AreEqual(twos, actual);
        }
        [TestMethod]
        public void UnaryMinusNegatesAllMatrixValues()
        {
            var actual = R.UnaryMinus();
            Assert.AreEqual(Z, actual.Add(R));
        }
        [TestMethod]
        public void ArrayLeftDivideWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                S = A.ArrayLeftDivide(S);
            });
        }
        [TestMethod]
        public void ArrayLeftDivideGivesCorrectResult()
        {
            var actual = A.ArrayLeftDivide(R);
            Assert.AreEqual(O, actual);
        }
        [TestMethod]
        public void ArrayLeftDivideEqualsWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                A.ArrayLeftDivideEquals(S);
            });
        }
        [TestMethod]
        public void ArrayLeftDivideEqualsGivesCorrectResult()
        {
            A.ArrayLeftDivideEquals(R);
            Assert.AreEqual(O, A);
        }
        [TestMethod]
        public void ArrayRightDivideWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                S = A.ArrayRightDivide(S);
            });
        }
        [TestMethod]
        public void ArrayRightDivideGivesCorrectResult()
        {
            var actual = A.ArrayRightDivide(R);
            Assert.AreEqual(O, actual);
        }
        [TestMethod]
        public void ArrayRightDivideEqualsWithNonConformantMatrixThrowsAnException()
        {
            Assert.ThrowsExactly<ArgumentException>(() =>
            {
                A.ArrayRightDivideEquals(S);
            });
        }
        [TestMethod]
        public void ArrayRightDivideEqualsGivesCorrectResult()
        {
            A.ArrayRightDivideEquals(R);
            Assert.AreEqual(O, A);
        }

    }
}
