﻿using System;
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
        GeneralMatrix A, B, C, Z, O, I, R, S, X, SUB, M, T, SQ, DEF, SOL;
        int errorCount = 0;
        int warningCount = 0;
        double tmp;
        double[] columnwise { get { return new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 }; } }
        double[] rowwise { get { return new[] { 1.0, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 3.0, 6.0, 9.0, 12.0 }; } } 
        double[][] avals { get{ return new[] { new[] { 1.0, 4.0, 7.0, 10.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } }; } }
 
        //double[][] rankdef = avals;
        double[][] tvals { get { return new[] {new[] {1.0, 2.0, 3.0}, new[] {4.0, 5.0, 6.0}, new[] {7.0, 8.0, 9.0}, new[] {10.0, 11.0, 12.0}}; }} 
        double[][] subavals { get { return new[] {new[] {5.0, 8.0, 11.0}, new[] {6.0, 9.0, 12.0}}; }} 
        double[][] rvals { get { return new[] {new[] {1.0, 4.0, 7.0}, new[] {2.0, 5.0, 8.0, 11.0}, new[] {3.0, 6.0, 9.0, 12.0}}; }}
        double[][] pvals { get { return new[] { new[] { 1.0, 1.0, 1.0 }, new[] { 1.0, 2.0, 3.0 }, new[] { 1.0, 3.0, 6.0 } }; } } 
        double[][] ivals { get { return new[] { new[] { 1.0, 0.0, 0.0, 0.0 }, new[] { 0.0, 1.0, 0.0, 0.0 }, new[] { 0.0, 0.0, 1.0, 0.0 } }; } } 
        double[][] evals = { new[] { 0.0, 1.0, 0.0, 0.0 }, new[] { 1.0, 0.0, 2e-7, 0.0 }, new[] { 0.0, -2e-7, 0.0, 1.0 }, new[] { 0.0, 0.0, 1.0, 0.0 } };
        double[][] square = { new[] { 166.0, 188.0, 210.0 }, new[] { 188.0, 214.0, 240.0 }, new[] { 210.0, 240.0, 270.0 } };
        double[][] sqSolution = { new[] { 13.0 }, new[] { 15.0 } };
        double[][] condmat = { new[] { 1.0, 3.0 }, new[] { 7.0, 9.0 } };
        int rows = 3, cols = 4;
        int invalidld = 5; /* should trigger bad shape for construction with val */
        int raggedr = 0; /* (raggedr,raggedc) should be out of bounds in ragged array */
        int raggedc = 4;
        int validld = 3; /* leading dimension of intended test Matrices */
        int nonconformld = 4; /* leading dimension which is valid, but nonconforming */
        int ib = 1, ie = 2, jb = 1, je = 3; /* index ranges for sub GeneralMatrix */
        int[] rowindexset = new int[] { 1, 2 };
        int[] badrowindexset = new int[] { 1, 3 };
        int[] columnindexset = new int[] { 1, 2, 3 };
        int[] badcolumnindexset = new int[] { 1, 2, 4 };
        double columnsummax = 33.0;
        double rowsummax = 30.0;
        double sumofdiagonals = 15;
        double sumofsquares = 650;

        /// <summary>Check norm of difference of Matrices. *</summary>

        private static void check(GeneralMatrix x, GeneralMatrix y)
        {
            double eps = Math.Pow(2.0, -52.0);
            if (x.Norm1() == 0.0 & y.Norm1() < 10 * eps)
                return;
            if (y.Norm1() == 0.0 & x.Norm1() < 10 * eps)
                return;
            if (x.Subtract(y).Norm1() > 1000 * eps * Math.Max(x.Norm1(), y.Norm1()))
            {
                throw new SystemException("The norm of (X-Y) is too large: " + x.Subtract(y).Norm1());
            }
        }
        /// <summary>
        /// check that exception is thrown in packed constructor with invalid length 
        /// </summary>
        [TestMethod][ExpectedException(typeof(ArgumentException))]
        public void ArgumentExceptionIsThrownInPackedConstructorWithInvalidLength()
        {
            int invalid = 5; /* should trigger bad shape for construction with val */
            var a = new GeneralMatrix(columnwise, invalid);
            Assert.Inconclusive(a.ToString());
        }
        /// <summary>
        /// check that exception is thrown in default constructor if input array is 'ragged' *
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void ArgumentExceptionIsThrownInConstructorWithRaggedInputArray()
        {
            double[][] array = { new[] { 1.0, 4.0, 7.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } };
            var a = new GeneralMatrix(array);
            Assert.Inconclusive(a.ToString());
        }
        /// <summary>
        /// check that exception is thrown in Create if input array is 'ragged' *
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void ArgumentExceptionIsThrownInCreateWithRaggedInputArray()
        {
            double[][] array = { new[] { 1.0, 4.0, 7.0 }, new[] { 2.0, 5.0, 8.0, 11.0 }, new[] { 3.0, 6.0, 9.0, 12.0 } };
            var a = GeneralMatrix.Create(array);
            Assert.Inconclusive(a.ToString());
        }

        [TestMethod]
        public void CanBeConstructedFromArrayAndRowOffset()
        {
            double[] packedColumns = new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            int numRows = 3; /* leading dimension of intended test Matrices */
            var a = new GeneralMatrix(packedColumns, numRows);
            Assert.IsTrue(a.RowDimension == numRows );
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
            Assert.AreEqual(3,matrix.RowDimension);
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
            Assert.AreEqual(expected, matrix.Array);
        }
        [TestMethod]
        public void ArrayCopyPropertyReturnsADeepCopyOfUnderlyingArray()
        {
            var matrix = new GeneralMatrix(avals);
            Assert.AreNotEqual(avals, matrix.ArrayCopy);
        }
        [TestMethod]
        public void ColumnPackedCopyPropertyReturnsAColumnPackedCopyOfUnderlyingArray()
        {
            var expected = new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
            var matrix = new GeneralMatrix(avals);
            Assert.IsTrue(expected.SequenceEqual(matrix.ColumnPackedCopy));
            matrix = new GeneralMatrix(expected,3);
            Assert.IsTrue(expected.SequenceEqual(matrix.ColumnPackedCopy));
            Assert.AreNotEqual(expected,matrix);
        }
        [TestMethod]
        public void RowPackedCopyPropertyReturnsARowPackedCopyOfUnderlyingArray()
        {
            var expected = new[] { 1.0, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 3.0, 6.0, 9.0, 12.0 };
            var matrix = new GeneralMatrix(avals);
            Assert.IsTrue(expected.SequenceEqual(matrix.RowPackedCopy));
            matrix = new GeneralMatrix(3, expected);
            Assert.IsTrue(expected.SequenceEqual(matrix.RowPackedCopy));
            Assert.AreNotEqual(expected, matrix);
        }
        [TestMethod]
        [ExpectedException(typeof(IndexOutOfRangeException))]
        public void ZeroElementAccessThrowsAnIndexOutOfRangeExceptionWhenColumnValueIsEqualToColumnDimension()
        {
            var matrix = new GeneralMatrix(avals);
            matrix.GetElement(matrix.RowDimension - 1, matrix.ColumnDimension);
        }
        [TestMethod]
        [ExpectedException(typeof(IndexOutOfRangeException))]
        public void ZeroElementAccessThrowsAnIndexOutOfRangeExceptionWhenRowValueIsEqualToRowDimension()
        {
            var matrix = new GeneralMatrix(avals);
            matrix.GetElement(matrix.RowDimension, matrix.ColumnDimension-1);
        }
        [TestMethod]
        public void ZeroBasedGetElementMethodReturnsCorrectElement()
        {
            var matrix = new GeneralMatrix(avals);
            var actual = matrix.GetElement(matrix.RowDimension-1, matrix.ColumnDimension-1);
            Assert.AreEqual(12.0,actual);
        }
        [TestMethod]
        [ExpectedException(typeof(IndexOutOfRangeException))]
        public void GetSubMatrixThrowAnIndexOutOfRangeExceptionWhenFinalRowIndexExceedsRowDimension()
        {
            var matrix = new GeneralMatrix(avals);
            matrix.GetMatrix(ib, ie + matrix.RowDimension + 1, jb, je);
        }
        [TestMethod]
        [ExpectedException(typeof(IndexOutOfRangeException))]
        public void GetSubMatrixThrowAnIndexOutOfRangeExceptionWhenFinalColumnIndexExceedsColumnDimension()
        {
            var matrix = new GeneralMatrix(avals);
            matrix.GetMatrix(ib, ie, jb, je + matrix.ColumnDimension + 1);
        }
        [TestMethod]
        public void GetSubMatrixGetsCorrectSubMatrix()
        {
            var matrix = new GeneralMatrix(avals);
            var expected = new GeneralMatrix(subavals);
            var actual = matrix.GetMatrix(ib, ie, jb, je);
            Assert.AreEqual(expected,actual);
        }
    }
}
