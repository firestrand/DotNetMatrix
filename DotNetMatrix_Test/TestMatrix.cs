using System;
using DotNetMatrix;
namespace DotNetMatrix_Test;

// SUMMARY:TestMatrix tests the functionality of the DotNetMatrix GeneralMatrix class and associated decompositions.
// <P>
// Run the test from the command line using
// <BLOCKQUOTE><PRE><CODE>
// DotNetMatrix.test.TestMatrix
// </CODE></PRE></BLOCKQUOTE>
// Detailed output is provided indicating the functionality being tested
// and whether the functionality is correctly implemented.   Exception handling
// is also tested.
// <P>
// The test is designed to run to completion and give a summary of any implementation errors
// encountered. The final output should be:
// <BLOCKQUOTE><PRE><CODE>
// TestMatrix completed.
// Total errors reported: n1
// Total warning reported: n2
// </CODE></PRE></BLOCKQUOTE>
// If the test does not run to completion, this indicates that there is a
// substantial problem within the implementation that was not anticipated in the test design.
// The stopping point should give an indication of where the problem exists.
//
//
public class TestMatrix
{
    [STAThread]
    public static void Main(string[] argv) => Run(argv);

    public static int Run(string[] argv)
    {
        GeneralMatrix A, B, C, Z, O, R, S, X, SUB, M, T, SQ, DEF, SOL;
        int errorCount = 0;
        int warningCount = 0;
        double[] columnwise = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
        double[] rowwise = new double[] { 1.0, 4.0, 7.0, 10.0, 2.0, 5.0, 8.0, 11.0, 3.0, 6.0, 9.0, 12.0 };
        double[][] avals = { new double[] { 1.0, 4.0, 7.0, 10.0 }, new double[] { 2.0, 5.0, 8.0, 11.0 }, new double[] { 3.0, 6.0, 9.0, 12.0 } };
        double[][] rankdef = avals;
        double[][] tvals = { new double[] { 1.0, 2.0, 3.0 }, new double[] { 4.0, 5.0, 6.0 }, new double[] { 7.0, 8.0, 9.0 }, new double[] { 10.0, 11.0, 12.0 } };
        double[][] subavals = { new double[] { 5.0, 8.0, 11.0 }, new double[] { 6.0, 9.0, 12.0 } };
        double[][] rvals = { new double[] { 1.0, 4.0, 7.0 }, new double[] { 2.0, 5.0, 8.0, 11.0 }, new double[] { 3.0, 6.0, 9.0, 12.0 } };
        double[][] pvals = { new double[] { 1.0, 1.0, 1.0 }, new double[] { 1.0, 2.0, 3.0 }, new double[] { 1.0, 3.0, 6.0 } };
        double[][] ivals = { new double[] { 1.0, 0.0, 0.0, 0.0 }, new double[] { 0.0, 1.0, 0.0, 0.0 }, new double[] { 0.0, 0.0, 1.0, 0.0 } };
        double[][] evals = { new double[] { 0.0, 1.0, 0.0, 0.0 }, new double[] { 1.0, 0.0, 2e-7, 0.0 }, new double[] { 0.0, -2e-7, 0.0, 1.0 }, new double[] { 0.0, 0.0, 1.0, 0.0 } };
        double[][] square = { new double[] { 166.0, 188.0, 210.0 }, new double[] { 188.0, 214.0, 240.0 }, new double[] { 210.0, 240.0, 270.0 } };
        double[][] sqSolution = { new double[] { 13.0 }, new double[] { 15.0 } };
        double[][] condmat = { new double[] { 1.0, 3.0 }, new double[] { 7.0, 9.0 } };
        int validld = 3; /* leading dimension of intended test Matrices */
        int nonconformld = 4; /* leading dimension which is valid, but nonconforming */
        int[] rowindexset = new int[] { 1, 2 };
        int[] badrowindexset = new int[] { 1, 3 };
        int[] columnindexset = new int[] { 1, 2, 3 };
        int[] badcolumnindexset = new int[] { 1, 2, 4 };
        double columnsummax = 33.0;
        double rowsummax = 30.0;
        double sumofdiagonals = 15;
        double sumofsquares = 650;

        // SUMMARY:Constructors and constructor-like methods:
        // double[], int
        // double[][]
        // int, int
        // int, int, double
        // int, int, double[][]
        // Create(double[][])
        // Random(int,int)
        // Identity(int)
        //
        //

        //Print("\nTesting constructors and constructor-like methods...\n");
        //try
        //{
        //    // SUMMARY:check that exception is thrown in packed constructor with invalid length *
        //    A = new GeneralMatrix(columnwise, invalidld);
        //    errorCount = TryFailure(errorCount, "Catch invalid length in packed constructor... ", "exception not thrown for invalid input");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("Catch invalid length in packed constructor... ", e.Message);
        //}
        //try
        //{
        //    // SUMMARY:check that exception is thrown in default constructor
        //    // if input array is 'ragged' *
        //    //
        //    A = new GeneralMatrix(rvals);
        //    tmp = A.GetElement(raggedr, raggedc);
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("Catch ragged input to default constructor... ", e.Message);
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "Catch ragged input to constructor... ", "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    // SUMMARY:check that exception is thrown in Create
        //    // if input array is 'ragged' *
        //    //
        //    A = GeneralMatrix.Create(rvals);
        //    tmp = A.GetElement(raggedr, raggedc);
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("Catch ragged input to Create... ", e.Message);
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "Catch ragged input to Create... ", "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
        //    System.Console.Out.WriteLine(e.Message);
        //}

        //A = new GeneralMatrix(columnwise, validld);
        //B = new GeneralMatrix(avals);
        //tmp = B.GetElement(0, 0);
        //avals[0][0] = 0.0;
        //C = B.Subtract(A);
        //avals[0][0] = tmp;
        //B = GeneralMatrix.Create(avals);
        //tmp = B.GetElement(0, 0);
        //avals[0][0] = 0.0;
        //if ((tmp - B.GetElement(0, 0)) != 0.0)
        //{
        //    // SUMMARY:check that Create behaves properly *
        //    errorCount = TryFailure(errorCount, "Create... ", "Copy not effected... data visible outside");
        //}
        //else
        //{
        //    TrySuccess("Create... ", "");
        //}
        //avals[0][0] = columnwise[0];

        //I = new GeneralMatrix(ivals);
        //try
        //{
        //    Check(I, GeneralMatrix.Identity(3, 4));
        //    TrySuccess("Identity... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "Identity... ", "Identity GeneralMatrix not successfully created");
        //    System.Console.Out.WriteLine(e.Message);
        //}

        // SUMMARY:Access Methods:
        // getColumnDimension()
        // getRowDimension()
        // getArray()
        // getArrayCopy()
        // getColumnPackedCopy()
        // getRowPackedCopy()
        // get(int,int)
        // GetMatrix(int,int,int,int)
        // GetMatrix(int,int,int[])
        // GetMatrix(int[],int,int)
        // GetMatrix(int[],int[])
        // set(int,int,double)
        // SetMatrix(int,int,int,int,GeneralMatrix)
        // SetMatrix(int,int,int[],GeneralMatrix)
        // SetMatrix(int[],int,int,GeneralMatrix)
        // SetMatrix(int[],int[],GeneralMatrix)
        //
        //

        //Print("\nTesting access methods...\n");

        // SUMMARY:Various get methods:
        //
        //

        B = new GeneralMatrix(avals);
        //if (B.RowDimension != rows)
        //{
        //    errorCount = TryFailure(errorCount, "getRowDimension... ", "");
        //}
        //else
        //{
        //    TrySuccess("getRowDimension... ", "");
        //}
        //if (B.ColumnDimension != cols)
        //{
        //    errorCount = TryFailure(errorCount, "getColumnDimension... ", "");
        //}
        //else
        //{
        //    TrySuccess("getColumnDimension... ", "");
        //}
        //B = new GeneralMatrix(avals);
        //double[][] barray = B.Array;
        //if (barray != avals)
        //{
        //    errorCount = TryFailure(errorCount, "getArray... ", "");
        //}
        //else
        //{
        //    TrySuccess("getArray... ", "");
        //}
        //barray = B.ArrayCopy;
        //if (barray == avals)
        //{
        //    errorCount = TryFailure(errorCount, "getArrayCopy... ", "data not (deep) copied");
        //}
        //try
        //{
        //    Check(barray, avals);
        //    TrySuccess("getArrayCopy... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "getArrayCopy... ", "data not successfully (deep) copied");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //double[] bpacked = B.ColumnPackedCopy;
        //try
        //{
        //    Check(bpacked, columnwise);
        //    TrySuccess("getColumnPackedCopy... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "getColumnPackedCopy... ", "data not successfully (deep) copied by columns");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //bpacked = B.RowPackedCopy;
        //try
        //{
        //    Check(bpacked, rowwise);
        //    TrySuccess("getRowPackedCopy... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "getRowPackedCopy... ", "data not successfully (deep) copied by rows");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    tmp = B.GetElement(B.RowDimension, B.ColumnDimension - 1);
        //    errorCount = TryFailure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        tmp = B.GetElement(B.RowDimension - 1, B.ColumnDimension);
        //        errorCount = TryFailure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("get(int,int)... OutofBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    if (B.GetElement(B.RowDimension - 1, B.ColumnDimension - 1) != avals[B.RowDimension - 1][B.ColumnDimension - 1])
        //    {
        //        errorCount = TryFailure(errorCount, "get(int,int)... ", "GeneralMatrix entry (i,j) not successfully retreived");
        //    }
        //    else
        //    {
        //        TrySuccess("get(int,int)... ", "");
        //    }
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "get(int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        SUB = new GeneralMatrix(subavals);
        //try
        //{
        //    M = B.GetMatrix(ib, ie + B.RowDimension + 1, jb, je);
        //    errorCount = TryFailure(errorCount, "GetMatrix(int,int,int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        M = B.GetMatrix(ib, ie, jb, je + B.ColumnDimension + 1);
        //        errorCount = TryFailure(errorCount, "GetMatrix(int,int,int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("GetMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int,int,int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    M = B.GetMatrix(ib, ie, jb, je);
        //    try
        //    {
        //        Check(SUB, M);
        //        TrySuccess("GetMatrix(int,int,int,int)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "GetMatrix(int,int,int,int)... ", "submatrix not successfully retreived");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int,int,int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e.Message);
        //}

        //try
        //{
        //    M = B.GetMatrix(ib, ie, badcolumnindexset);
        //    errorCount = TryFailure(errorCount, "GetMatrix(int,int,int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        M = B.GetMatrix(ib, ie + B.RowDimension + 1, columnindexset);
        //        errorCount = TryFailure(errorCount, "GetMatrix(int,int,int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("GetMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int,int,int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    M = B.GetMatrix(ib, ie, columnindexset);
        //    try
        //    {
        //        Check(SUB, M);
        //        TrySuccess("GetMatrix(int,int,int[])... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "GetMatrix(int,int,int[])... ", "submatrix not successfully retreived");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int,int,int[])... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    M = B.GetMatrix(badrowindexset, jb, je);
        //    errorCount = TryFailure(errorCount, "GetMatrix(int[],int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        M = B.GetMatrix(rowindexset, jb, je + B.ColumnDimension + 1);
        //        errorCount = TryFailure(errorCount, "GetMatrix(int[],int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("GetMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int[],int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    M = B.GetMatrix(rowindexset, jb, je);
        //    try
        //    {
        //        Check(SUB, M);
        //        TrySuccess("GetMatrix(int[],int,int)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "GetMatrix(int[],int,int)... ", "submatrix not successfully retreived");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int[],int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    M = B.GetMatrix(badrowindexset, columnindexset);
        //    errorCount = TryFailure(errorCount, "GetMatrix(int[],int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        M = B.GetMatrix(rowindexset, badcolumnindexset);
        //        errorCount = TryFailure(errorCount, "GetMatrix(int[],int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("GetMatrix(int[],int[])... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int[],int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    M = B.GetMatrix(rowindexset, columnindexset);
        //    try
        //    {
        //        Check(SUB, M);
        //        TrySuccess("GetMatrix(int[],int[])... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "GetMatrix(int[],int[])... ", "submatrix not successfully retreived");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    errorCount = TryFailure(errorCount, "GetMatrix(int[],int[])... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e.Message);
        //}

        // SUMMARY:Various set methods:
        //
        //

        //try
        //{
        //    B.SetElement(B.RowDimension, B.ColumnDimension - 1, 0.0);
        //    errorCount = TryFailure(errorCount, "set(int,int,double)... ", "OutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        B.SetElement(B.RowDimension - 1, B.ColumnDimension, 0.0);
        //        errorCount = TryFailure(errorCount, "set(int,int,double)... ", "OutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("set(int,int,double)... OutofBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "set(int,int,double)... ", "OutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetElement(ib, jb, 0.0);
        //    tmp = B.GetElement(ib, jb);
        //    try
        //    {
        //        Check(tmp, 0.0);
        //        TrySuccess("set(int,int,double)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "set(int,int,double)... ", "GeneralMatrix element not successfully set");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //}
        //catch (System.IndexOutOfRangeException e1)
        //{
        //    errorCount = TryFailure(errorCount, "set(int,int,double)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        M = new GeneralMatrix(2, 3, 0.0);
        //try
        //{
        //    B.SetMatrix(ib, ie + B.RowDimension + 1, jb, je, M);
        //    errorCount = TryFailure(errorCount, "SetMatrix(int,int,int,int,GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        B.SetMatrix(ib, ie, jb, je + B.ColumnDimension + 1, M);
        //        errorCount = TryFailure(errorCount, "SetMatrix(int,int,int,int,GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("SetMatrix(int,int,int,int,GeneralMatrix)... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int,int,int,int,GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(ib, ie, jb, je, M);
        //    try
        //    {
        //        Check(M.Subtract(B.GetMatrix(ib, ie, jb, je)), M);
        //        TrySuccess("SetMatrix(int,int,int,int,GeneralMatrix)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "SetMatrix(int,int,int,int,GeneralMatrix)... ", "submatrix not successfully set");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //    B.SetMatrix(ib, ie, jb, je, SUB);
        //}
        //catch (System.IndexOutOfRangeException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int,int,int,int,GeneralMatrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(ib, ie + B.RowDimension + 1, columnindexset, M);
        //    errorCount = TryFailure(errorCount, "SetMatrix(int,int,int[],GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        B.SetMatrix(ib, ie, badcolumnindexset, M);
        //        errorCount = TryFailure(errorCount, "SetMatrix(int,int,int[],GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("SetMatrix(int,int,int[],GeneralMatrix)... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int,int,int[],GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(ib, ie, columnindexset, M);
        //    try
        //    {
        //        Check(M.Subtract(B.GetMatrix(ib, ie, columnindexset)), M);
        //        TrySuccess("SetMatrix(int,int,int[],GeneralMatrix)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "SetMatrix(int,int,int[],GeneralMatrix)... ", "submatrix not successfully set");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //    B.SetMatrix(ib, ie, jb, je, SUB);
        //}
        //catch (System.IndexOutOfRangeException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int,int,int[],GeneralMatrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(rowindexset, jb, je + B.ColumnDimension + 1, M);
        //    errorCount = TryFailure(errorCount, "SetMatrix(int[],int,int,GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        B.SetMatrix(badrowindexset, jb, je, M);
        //        errorCount = TryFailure(errorCount, "SetMatrix(int[],int,int,GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("SetMatrix(int[],int,int,GeneralMatrix)... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int[],int,int,GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(rowindexset, jb, je, M);
        //    try
        //    {
        //        Check(M.Subtract(B.GetMatrix(rowindexset, jb, je)), M);
        //        TrySuccess("SetMatrix(int[],int,int,GeneralMatrix)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "SetMatrix(int[],int,int,GeneralMatrix)... ", "submatrix not successfully set");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //    B.SetMatrix(ib, ie, jb, je, SUB);
        //}
        //catch (System.IndexOutOfRangeException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int[],int,int,GeneralMatrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(rowindexset, badcolumnindexset, M);
        //    errorCount = TryFailure(errorCount, "SetMatrix(int[],int[],GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //}
        //catch (System.IndexOutOfRangeException e)
        //{
        //    System.Console.Out.WriteLine(e.Message);
        //    try
        //    {
        //        B.SetMatrix(badrowindexset, columnindexset, M);
        //        errorCount = TryFailure(errorCount, "SetMatrix(int[],int[],GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    }
        //    catch (System.IndexOutOfRangeException e1)
        //    {
        //        TrySuccess("SetMatrix(int[],int[],GeneralMatrix)... ArrayIndexOutOfBoundsException... ", "");
        //        System.Console.Out.WriteLine(e1.Message);
        //    }
        //}
        //catch (System.ArgumentException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int[],int[],GeneralMatrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
        //    System.Console.Out.WriteLine(e1.Message);
        //}
        //try
        //{
        //    B.SetMatrix(rowindexset, columnindexset, M);
        //    try
        //    {
        //        Check(M.Subtract(B.GetMatrix(rowindexset, columnindexset)), M);
        //        TrySuccess("SetMatrix(int[],int[],GeneralMatrix)... ", "");
        //    }
        //    catch (System.SystemException e)
        //    {
        //        errorCount = TryFailure(errorCount, "SetMatrix(int[],int[],GeneralMatrix)... ", "submatrix not successfully set");
        //        System.Console.Out.WriteLine(e.Message);
        //    }
        //}
        //catch (System.IndexOutOfRangeException e1)
        //{
        //    errorCount = TryFailure(errorCount, "SetMatrix(int[],int[],GeneralMatrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
        //    System.Console.Out.WriteLine(e1.Message);
        //}

        // SUMMARY:Array-like methods:
        // Subtract
        // SubtractEquals
        // Add
        // AddEquals
        // ArrayLeftDivide
        // ArrayLeftDivideEquals
        // ArrayRightDivide
        // ArrayRightDivideEquals
        // arrayTimes
        // ArrayMultiplyEquals
        // uminus
        //
        //
        A = new GeneralMatrix(columnwise, validld);
        Print("\nTesting array-like methods...\n");
        S = new GeneralMatrix(columnwise, nonconformld);
        R = GeneralMatrix.Random(A.RowDimension, A.ColumnDimension);
        A = R;
        //try
        //{
        //    S = A.Subtract(S);
        //    errorCount = TryFailure(errorCount, "Subtract conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("Subtract conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //if (A.Subtract(R).Norm1() != 0.0)
        //{
        //    errorCount = TryFailure(errorCount, "Subtract... ", "(difference of identical Matrices is nonzero,\nSubsequent use of Subtract should be suspect)");
        //}
        //else
        //{
        //    TrySuccess("Subtract... ", "");
        //}
        A = R.Copy();
        A.SubtractEquals(R);
        Z = new GeneralMatrix(A.RowDimension, A.ColumnDimension);
        //try
        //{
        //    A.SubtractEquals(S);
        //    errorCount = TryFailure(errorCount, "SubtractEquals conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("SubtractEquals conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //if (A.Subtract(Z).Norm1() != 0.0)
        //{
        //    errorCount = TryFailure(errorCount, "SubtractEquals... ", "(difference of identical Matrices is nonzero,\nSubsequent use of Subtract should be suspect)");
        //}
        //else
        //{
        //    TrySuccess("SubtractEquals... ", "");
        //}

        A = R.Copy();
        B = GeneralMatrix.Random(A.RowDimension, A.ColumnDimension);
        C = A.Subtract(B);
        //try
        //{
        //    S = A.Add(S);
        //    errorCount = TryFailure(errorCount, "Add conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("Add conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    Check(C.Add(B), A);
        //    TrySuccess("Add... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "Add... ", "(C = A - B, but C + B != A)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        C = A.Subtract(B);
        C.AddEquals(B);
        //try
        //{
        //    A.AddEquals(S);
        //    errorCount = TryFailure(errorCount, "AddEquals conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("AddEquals conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    Check(C, A);
        //    TrySuccess("AddEquals... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "AddEquals... ", "(C = A - B, but C = C + B != A)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //A = R.UnaryMinus();
        //try
        //{
        //    Check(A.Add(R), Z);
        //    TrySuccess("UnaryMinus... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "uminus... ", "(-A + A != zeros)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        A = R.Copy();
        O = new GeneralMatrix(A.RowDimension, A.ColumnDimension, 1.0);
        C = A.ArrayLeftDivide(R);
        //try
        //{
        //    S = A.ArrayLeftDivide(S);
        //    errorCount = TryFailure(errorCount, "ArrayLeftDivide conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("ArrayLeftDivide conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    Check(C, O);
        //    TrySuccess("ArrayLeftDivide... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "ArrayLeftDivide... ", "(M.\\M != ones)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    A.ArrayLeftDivideEquals(S);
        //    errorCount = TryFailure(errorCount, "ArrayLeftDivideEquals conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("ArrayLeftDivideEquals conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //A.ArrayLeftDivideEquals(R);
        //try
        //{
        //    Check(A, O);
        //    TrySuccess("ArrayLeftDivideEquals... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "ArrayLeftDivideEquals... ", "(M.\\M != ones)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //A = R.Copy();
        //try
        //{
        //    A.ArrayRightDivide(S);
        //    errorCount = TryFailure(errorCount, "ArrayRightDivide conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("ArrayRightDivide conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //C = A.ArrayRightDivide(R);
        //try
        //{
        //    Check(C, O);
        //    TrySuccess("ArrayRightDivide... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "ArrayRightDivide... ", "(M./M != ones)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //try
        //{
        //    A.ArrayRightDivideEquals(S);
        //    errorCount = TryFailure(errorCount, "ArrayRightDivideEquals conformance check... ", "nonconformance not raised");
        //}
        //catch (System.ArgumentException e)
        //{
        //    TrySuccess("ArrayRightDivideEquals conformance check... ", "");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        //A.ArrayRightDivideEquals(R);
        //try
        //{
        //    Check(A, O);
        //    TrySuccess("ArrayRightDivideEquals... ", "");
        //}
        //catch (System.SystemException e)
        //{
        //    errorCount = TryFailure(errorCount, "ArrayRightDivideEquals... ", "(M./M != ones)");
        //    System.Console.Out.WriteLine(e.Message);
        //}
        A = R.Copy();
        B = GeneralMatrix.Random(A.RowDimension, A.ColumnDimension);
        try
        {
            S = A.ArrayMultiply(S);
            errorCount = TryFailure(errorCount, "arrayTimes conformance check... ", "nonconformance not raised");
        }
        catch (System.ArgumentException e)
        {
            TrySuccess("arrayTimes conformance check... ", "");
            System.Console.Out.WriteLine(e.Message);
        }
        C = A.ArrayMultiply(B);
        try
        {
            Check(C.ArrayRightDivideEquals(B), A);
            TrySuccess("arrayTimes... ", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "arrayTimes... ", "(A = R, C = A.*B, but C./B != A)");
            System.Console.Out.WriteLine(e.Message);
        }
        try
        {
            A.ArrayMultiplyEquals(S);
            errorCount = TryFailure(errorCount, "ArrayMultiplyEquals conformance check... ", "nonconformance not raised");
        }
        catch (System.ArgumentException e)
        {
            TrySuccess("ArrayMultiplyEquals conformance check... ", "");
            System.Console.Out.WriteLine(e.Message);
        }
        A.ArrayMultiplyEquals(B);
        try
        {
            Check(A.ArrayRightDivideEquals(B), R);
            TrySuccess("ArrayMultiplyEquals... ", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "ArrayMultiplyEquals... ", "(A = R, A = A.*B, but A./B != R)");
            System.Console.Out.WriteLine(e.Message);
        }

        // SUMMARY:LA methods:
        // Transpose
        // Multiply
        // Condition
        // Rank
        // Determinant
        // trace
        // Norm1
        // norm2
        // normF
        // normInf
        // Solve
        // solveTranspose
        // Inverse
        // chol
        // Eigen
        // lu
        // qr
        // svd
        //
        //

        Print("\nTesting linear algebra methods...\n");
        A = new GeneralMatrix(columnwise, 3);
        T = new GeneralMatrix(tvals);
        T = A.Transpose();
        try
        {
            Check(A.Transpose(), T);
            TrySuccess("Transpose...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Transpose()...", "Transpose unsuccessful");
            System.Console.Out.WriteLine(e.Message);
        }
        A.Transpose();
        try
        {
            Check(A.Norm1(), columnsummax);
            TrySuccess("Norm1...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Norm1()...", "incorrect norm calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        try
        {
            Check(A.NormInf(), rowsummax);
            TrySuccess("normInf()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "normInf()...", "incorrect norm calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        try
        {
            Check(A.NormF(), System.Math.Sqrt(sumofsquares));
            TrySuccess("normF...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "normF()...", "incorrect norm calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        try
        {
            Check(A.Trace(), sumofdiagonals);
            TrySuccess("trace()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "trace()...", "incorrect trace calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        try
        {
            Check(A.GetMatrix(0, A.RowDimension - 1, 0, A.RowDimension - 1).Determinant(), 0.0);
            TrySuccess("Determinant()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Determinant()...", "incorrect determinant calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        SQ = new GeneralMatrix(square);
        try
        {
            Check(A.Multiply(A.Transpose()), SQ);
            TrySuccess("Multiply(GeneralMatrix)...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Multiply(GeneralMatrix)...", "incorrect GeneralMatrix-GeneralMatrix product calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        try
        {
            Check(A.Multiply(0.0), Z);
            TrySuccess("Multiply(double)...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Multiply(double)...", "incorrect GeneralMatrix-scalar product calculation");
            System.Console.Out.WriteLine(e.Message);
        }

        A = new GeneralMatrix(columnwise, 4);
        QRDecomposition QR = A.Qrd();
        R = QR.R;
        try
        {
            Check(A, QR.Q.Multiply(R));
            TrySuccess("QRDecomposition...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "QRDecomposition...", "incorrect QR decomposition calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        SingularValueDecomposition SVD = A.Svd();
        try
        {
            Check(A, SVD.GetU().Multiply(SVD.S.Multiply(SVD.GetV().Transpose())));
            TrySuccess("SingularValueDecomposition...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "SingularValueDecomposition...", "incorrect singular value decomposition calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        DEF = new GeneralMatrix(rankdef);
        try
        {
            Check(DEF.Rank(), System.Math.Min(DEF.RowDimension, DEF.ColumnDimension) - 1);
            TrySuccess("Rank()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Rank()...", "incorrect Rank calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        B = new GeneralMatrix(condmat);
        SVD = B.Svd();
        double[] singularvalues = SVD.SingularValues;
        try
        {
            Check(B.Condition(), singularvalues[0] / singularvalues[System.Math.Min(B.RowDimension, B.ColumnDimension) - 1]);
            TrySuccess("Condition()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Condition()...", "incorrect condition number calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        int n = A.ColumnDimension;
        A = A.GetMatrix(0, n - 1, 0, n - 1);
        A.SetElement(0, 0, 0.0);
        LUDecomposition LU = A.Lud();
        try
        {
            Check(A.GetMatrix(LU.Pivot, 0, n - 1), LU.L.Multiply(LU.U));
            TrySuccess("LUDecomposition...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "LUDecomposition...", "incorrect LU decomposition calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        X = A.Inverse();
        try
        {
            Check(A.Multiply(X), GeneralMatrix.Identity(3, 3));
            TrySuccess("Inverse()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Inverse()...", "incorrect Inverse calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        O = new GeneralMatrix(SUB.RowDimension, 1, 1.0);
        SOL = new GeneralMatrix(sqSolution);
        SQ = SUB.GetMatrix(0, SUB.RowDimension - 1, 0, SUB.RowDimension - 1);
        try
        {
            Check(SQ.Solve(SOL), O);
            TrySuccess("Solve()...", "");
        }
        catch (System.ArgumentException e1)
        {
            errorCount = TryFailure(errorCount, "Solve()...", e1.Message);
            System.Console.Out.WriteLine(e1.Message);
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "Solve()...", e.Message);
            System.Console.Out.WriteLine(e.Message);
        }
        A = new GeneralMatrix(pvals);
        CholeskyDecomposition Chol = A.Chol();
        GeneralMatrix L = Chol.GetL();
        try
        {
            Check(A, L.Multiply(L.Transpose()));
            TrySuccess("CholeskyDecomposition...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "CholeskyDecomposition...", "incorrect Cholesky decomposition calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        X = Chol.Solve(GeneralMatrix.Identity(3, 3));
        try
        {
            Check(A.Multiply(X), GeneralMatrix.Identity(3, 3));
            TrySuccess("CholeskyDecomposition Solve()...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "CholeskyDecomposition Solve()...", "incorrect Choleskydecomposition Solve calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        EigenvalueDecomposition Eig = A.Eigen();
        GeneralMatrix D = Eig.D;
        GeneralMatrix V = Eig.GetV();
        try
        {
            Check(A.Multiply(V), V.Multiply(D));
            TrySuccess("EigenvalueDecomposition (symmetric)...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "EigenvalueDecomposition (symmetric)...", "incorrect symmetric Eigenvalue decomposition calculation");
            System.Console.Out.WriteLine(e.Message);
        }
        A = new GeneralMatrix(evals);
        Eig = A.Eigen();
        D = Eig.D;
        V = Eig.GetV();
        try
        {
            Check(A.Multiply(V), V.Multiply(D));
            TrySuccess("EigenvalueDecomposition (nonsymmetric)...", "");
        }
        catch (System.SystemException e)
        {
            errorCount = TryFailure(errorCount, "EigenvalueDecomposition (nonsymmetric)...", "incorrect nonsymmetric Eigenvalue decomposition calculation");
            System.Console.Out.WriteLine(e.Message);
        }

        Print("\nTestMatrix completed.\n");
        Print("Total errors reported: " + System.Convert.ToString(errorCount) + "\n");
        Print("Total warnings reported: " + System.Convert.ToString(warningCount) + "\n");
        return errorCount;
    }

    // SUMMARY:private utility routines *

    // SUMMARY:Check magnitude of difference of scalars. *

    private static void Check(double x, double y)
    {
        double eps = System.Math.Pow(2.0, -52.0);
        if (x == 0 & System.Math.Abs(y) < 10 * eps)
            return;
        if (y == 0 & System.Math.Abs(x) < 10 * eps)
            return;
        if (System.Math.Abs(x - y) > 10 * eps * System.Math.Max(System.Math.Abs(x), System.Math.Abs(y)))
        {
            throw new System.SystemException("The difference x-y is too large: x = " + x.ToString() + "  y = " + y.ToString());
        }
    }

    // SUMMARY:Check norm of difference of "vectors". *

    internal static void Check(double[] x, double[] y)
    {
        if (x.Length == y.Length)
        {
            for (int i = 0; i < x.Length; i++)
            {
                Check(x[i], y[i]);
            }
        }
        else
        {
            throw new System.SystemException("Attempt to compare vectors of different lengths");
        }
    }

    // SUMMARY:Check norm of difference of arrays. *

    internal static void Check(double[][] x, double[][] y)
    {
        GeneralMatrix A = new GeneralMatrix(x);
        GeneralMatrix B = new GeneralMatrix(y);
        Check(A, B);
    }

    // SUMMARY:Check norm of difference of Matrices. *

    private static void Check(GeneralMatrix X, GeneralMatrix Y)
    {
        double eps = System.Math.Pow(2.0, -52.0);
        if (X.Norm1() == 0.0 & Y.Norm1() < 10 * eps)
            return;
        if (Y.Norm1() == 0.0 & X.Norm1() < 10 * eps)
            return;
        if (X.Subtract(Y).Norm1() > 1000 * eps * System.Math.Max(X.Norm1(), Y.Norm1()))
        {
            throw new System.SystemException("The norm of (X-Y) is too large: " + X.Subtract(Y).Norm1().ToString());
        }
    }

    // SUMMARY:Shorten spelling of print. *

    private static void Print(System.String s)
    {
        System.Console.Out.Write(s);
    }

    // SUMMARY:Print appropriate messages for successful outcome try *

    private static void TrySuccess(System.String s, System.String e)
    {
        Print(">    " + s + "success\n");
        if ((System.Object)e != (System.Object)"")
        {
            Print(">      Message: " + e + "\n");
        }
    }
    // SUMMARY:Print appropriate messages for unsuccessful outcome try *

    private static int TryFailure(int count, System.String s, System.String e)
    {
        Print(">    " + s + "*** failure ***\n>      Message: " + e + "\n");
        return ++count;
    }

    // SUMMARY:Print appropriate messages for unsuccessful outcome try *

    internal static int TryWarning(int count, System.String s, System.String e)
    {
        Print(">    " + s + "*** warning ***\n>      Message: " + e + "\n");
        return ++count;
    }
}
