using System;
using System.Runtime.Serialization;

namespace DotNetMatrix;

/// <summary>
/// A mutable dense matrix of real double-precision values with arithmetic, submatrix,
/// linear-system, and matrix-decomposition operations.
/// </summary>
/// <remarks>
/// Cholesky and eigenvalue decompositions require square matrices. LU and QR support square
/// and tall matrices, while singular value decomposition also supports wide matrices.
/// The <see cref="Array"/> property exposes borrowed mutable storage; use <see cref="ArrayCopy"/>
/// or <see cref="Copy"/> for an independent snapshot. Callers must synchronize concurrent mutations.
/// </remarks>
/// <example>
/// Solve A*x = b and compute the infinity norm of the residual:
/// <code>
/// double[][] values = [[1, 2, 3], [4, 5, 6], [7, 8, 10]];
/// GeneralMatrix a = new(values);
/// GeneralMatrix b = GeneralMatrix.Random(3, 1);
/// GeneralMatrix x = a.Solve(b);
/// double residualNorm = a.Multiply(x).Subtract(b).NormInf();
/// </code>
/// </example>
/// <author>
///   The MathWorks, Inc. and the National Institute of Standards and Technology.
/// </author>
/// <version>  5 August 1998
/// </version>
[Serializable]
public partial class GeneralMatrix : ICloneable, ISerializable, IDisposable, IEquatable<GeneralMatrix>
{
    #region Class variables

    /// <summary>
    ///   Array for internal storage of elements.
    ///   @serial internal array storage.
    /// </summary>
    private readonly double[][] _a;

    /// <summary>
    ///   Row and column dimensions.
    ///   @serial row dimension.
    ///   @serial column dimension.
    /// </summary>
    private readonly int _m;

    /// <summary>
    ///   Row and column dimensions.
    ///   @serial row dimension.
    ///   @serial column dimension.
    /// </summary>
    private readonly int _n;

    #endregion //  Class variables

    #region Constructors

    /// <summary>
    ///   Construct an m-by-n matrix of zeros.
    /// </summary>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    public GeneralMatrix(int m, int n)
    {
        _m = m;
        _n = n;
        _a = new double[m][];
        for (int i = 0; i < m; i++)
        {
            _a[i] = new double[n];
        }
    }

    /// <summary>
    ///   Construct an m-by-n constant matrix.
    /// </summary>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    /// <param name = "s">   Fill the matrix with this scalar value.
    /// </param>
    public GeneralMatrix(int m, int n, double s)
    {
        _m = m;
        _n = n;
        _a = new double[m][];
        for (int i = 0; i < m; i++)
        {
            _a[i] = new double[n];
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                _a[i][j] = s;
            }
        }
    }

    /// <summary>
    ///   Construct a matrix from a 2-D array.
    /// </summary>
    /// <param name = "a">   Two-dimensional array of doubles.
    /// </param>
    /// <exception cref = "System.ArgumentException">   All rows must have the same length
    /// </exception>
    /// <seealso cref = "Create">
    /// </seealso>
    public GeneralMatrix(double[][] a)
    {
        _m = a.Length;
        _n = a[0].Length;
        for (int i = 0; i < _m; i++)
        {
            if (a[i].Length != _n)
            {
                throw new ArgumentException("All rows must have the same length.");
            }
        }
        _a = a;
    }

    /// <summary>
    ///   Construct a matrix quickly without checking arguments.
    /// </summary>
    /// <param name = "a">   Two-dimensional array of doubles.
    /// </param>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    public GeneralMatrix(double[][] a, int m, int n)
    {
        _a = a;
        _m = m;
        _n = n;
    }

    /// <summary>
    ///   Construct a matrix from a one-dimensional packed array
    /// </summary>
    /// <param name = "vals">One-dimensional array of doubles, packed by columns (ala Fortran).
    /// </param>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <exception cref = "System.ArgumentException">   Array length must be a multiple of m.
    /// </exception>
    public GeneralMatrix(double[] vals, int m)
    {
        _m = m;
        _n = (m != 0 ? vals.Length / m : 0);
        if (m * _n != vals.Length)
        {
            throw new ArgumentException("Array length must be a multiple of m.");
        }
        _a = new double[m][];
        for (int i = 0; i < m; i++)
        {
            _a[i] = new double[_n];
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = vals[i + j * m];
            }
        }
    }
    /// <summary>
    /// Copies a column-major block from a packed array into a new matrix.
    /// </summary>
    /// <param name="vals">The packed source values; the matrix does not retain the array.</param>
    /// <param name="m">The number of rows in the result.</param>
    /// <param name="n">The number of columns in the result.</param>
    /// <param name="startIndex">The first element of the column-major block in <paramref name="vals"/>.</param>
    /// <exception cref="ArgumentException">The remaining source does not contain the requested block.</exception>
    public GeneralMatrix(double[] vals, int m, int n, int startIndex)
    {
        _m = m;
        _n = n;
        if (!(_m * _n <= vals.Length - startIndex))
        {
            throw new ArgumentException("Array length must be a multiple of m and sized correctly.");
        }
        double[] subArray = new double[_m * _n];
        System.Array.Copy(vals, startIndex, subArray, 0, _m * _n);

        _a = new double[_m][];
        for (int i = 0; i < _m; i++)
        {
            _a[i] = new double[_n];
        }
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = subArray[i + j * _m];
            }
        }

    }
    /// <summary>
    ///   Construct a matrix from a one-dimensional packed array
    /// </summary>
    /// <param name = "vals">One-dimensional array of doubles, packed by rows.
    /// </param>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <exception cref = "System.ArgumentException">   Array length must be a multiple of m.
    /// </exception>
    public GeneralMatrix(int m, double[] vals)
    {
        _m = m;
        _n = (m != 0 ? vals.Length / m : 0);
        if (m * _n != vals.Length)
        {
            throw new ArgumentException("Array length must be a multiple of m.");
        }
        _a = new double[m][];
        for (int i = 0; i < m; i++)
        {
            _a[i] = new double[_n];
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = vals[i * _n + j];
            }
        }
    }
    #endregion //  Constructors

    #region Public Properties

    /// <summary>
    ///   Access the internal two-dimensional array.
    /// </summary>
    /// <remarks>Mutations through this borrowed array change the matrix. Prefer <see cref="ArrayCopy"/> when isolation is required.</remarks>
    public virtual double[][] Array
    {
        get { return _a; }
    }

    /// <summary>
    ///   Copy the internal two-dimensional array.
    /// </summary>
    /// <returns>     Two-dimensional array copy of matrix elements.
    /// </returns>
    public virtual double[][] ArrayCopy
    {
        get
        {
            var c = new double[_m][];
            for (int i = 0; i < _m; i++)
            {
                c[i] = new double[_n];
            }
            for (int i = 0; i < _m; i++)
            {
                for (int j = 0; j < _n; j++)
                {
                    c[i][j] = _a[i][j];
                }
            }
            return c;
        }
    }

    /// <summary>
    ///   Make a one-dimensional column packed copy of the internal array.
    /// </summary>
    /// <returns>     Matrix elements packed in a one-dimensional array by columns.
    /// </returns>
    public virtual double[] ColumnPackedCopy
    {
        get
        {
            var vals = new double[_m * _n];
            for (int i = 0; i < _m; i++)
            {
                for (int j = 0; j < _n; j++)
                {
                    vals[i + j * _m] = _a[i][j];
                }
            }
            return vals;
        }
    }

    /// <summary>
    ///   Make a one-dimensional row packed copy of the internal array.
    /// </summary>
    /// <returns>     Matrix elements packed in a one-dimensional array by rows.
    /// </returns>
    public virtual double[] RowPackedCopy
    {
        get
        {
            var vals = new double[_m * _n];
            for (int i = 0; i < _m; i++)
            {
                for (int j = 0; j < _n; j++)
                {
                    vals[i * _n + j] = _a[i][j];
                }
            }
            return vals;
        }
    }

    /// <summary>
    ///   Get row dimension.
    /// </summary>
    /// <returns>     m, the number of rows.
    /// </returns>
    public virtual int RowDimension
    {
        get { return _m; }
    }

    /// <summary>
    ///   Get column dimension.
    /// </summary>
    /// <returns>     n, the number of columns.
    /// </returns>
    public virtual int ColumnDimension
    {
        get { return _n; }
    }

    #endregion   // Public Properties

    #region	 Public Methods

    /// <summary>
    ///   Construct a matrix from a copy of a 2-D array.
    /// </summary>
    /// <param name = "a">   Two-dimensional array of doubles.
    /// </param>
    /// <exception cref = "System.ArgumentException">   All rows must have the same length
    /// </exception>
    public static GeneralMatrix Create(double[][] a)
    {
        int m = a.Length;
        int n = a[0].Length;
        var x = new GeneralMatrix(m, n);
        double[][] c = x.Array;
        for (int i = 0; i < m; i++)
        {
            if (a[i].Length != n)
            {
                throw new ArgumentException("All rows must have the same length.");
            }
            for (int j = 0; j < n; j++)
            {
                c[i][j] = a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   Make a deep copy of a matrix
    /// </summary>
    public virtual GeneralMatrix Copy()
    {
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        // Small rows do not amortize the bulk-copy call; retain the scalar path.
        if (_n < 16)
        {
            for (int i = 0; i < _m; i++)
            {
                for (int j = 0; j < _n; j++)
                {
                    c[i][j] = _a[i][j];
                }
            }
        }
        else
        {
            for (int i = 0; i < _m; i++)
            {
                double[] sourceRow = _a[i];
                double[] destinationRow = c[i];
                if (sourceRow.Length >= _n)
                {
                    System.Array.Copy(sourceRow, destinationRow, _n);
                }
                else
                {
                    // Preserve the legacy exception for malformed borrowed rows.
                    for (int j = 0; j < _n; j++)
                    {
                        destinationRow[j] = sourceRow[j];
                    }
                }
            }
        }
        return x;
    }

    /// <summary>
    ///   Get a single element.
    /// </summary>
    /// <param name = "i">   Row index.
    /// </param>
    /// <param name = "j">   Column index.
    /// </param>
    /// <returns>     A(i,j)
    /// </returns>
    /// <exception cref = "System.IndexOutOfRangeException">
    /// </exception>
    public virtual double GetElement(int i, int j)
    {
        return _a[i][j];
    }

    /// <summary>
    ///   Get a submatrix.
    /// </summary>
    /// <param name = "startRowIndex">  Initial row index
    /// </param>
    /// <param name = "endRowIndex">  Final row index
    /// </param>
    /// <param name = "j0">  Initial column index
    /// </param>
    /// <param name = "j1">  Final column index
    /// </param>
    /// <returns>     A(startRowIndex:endRowIndex,j0:j1)
    /// </returns>
    /// <exception cref = "System.IndexOutOfRangeException">   Submatrix indices
    /// </exception>
    public virtual GeneralMatrix GetMatrix(int startRowIndex, int endRowIndex, int j0, int j1)
    {
        var x = new GeneralMatrix(endRowIndex - startRowIndex + 1, j1 - j0 + 1);
        double[][] b = x.Array;
        try
        {
            for (int i = startRowIndex; i <= endRowIndex; i++)
            {
                for (int j = j0; j <= j1; j++)
                {
                    b[i - startRowIndex][j - j0] = _a[i][j];
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
        return x;
    }

    /// <summary>
    ///   Get a submatrix.
    /// </summary>
    /// <param name = "r">   Array of row indices.
    /// </param>
    /// <param name = "c">   Array of column indices.
    /// </param>
    /// <returns>     A(r(:),columnIndexes(:))
    /// </returns>
    /// <exception cref = "System.IndexOutOfRangeException">   Submatrix indices
    /// </exception>
    public virtual GeneralMatrix GetMatrix(int[] r, int[] c)
    {
        var x = new GeneralMatrix(r.Length, c.Length);
        double[][] b = x.Array;
        try
        {
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = 0; j < c.Length; j++)
                {
                    b[i][j] = _a[r[i]][c[j]];
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
        return x;
    }

    /// <summary>
    ///   Get a submatrix.
    /// </summary>
    /// <param name = "startRowIndex">  Initial row index
    /// </param>
    /// <param name = "endRowIndex">  Final row index
    /// </param>
    /// <param name = "columnIndexes">   Array of column indices.
    /// </param>
    /// <returns>     A(startRowIndex:endRowIndex,columnIndexes(:))
    /// </returns>
    /// <exception cref = "System.IndexOutOfRangeException">   Submatrix indices
    /// </exception>
    public virtual GeneralMatrix GetMatrix(int startRowIndex, int endRowIndex, int[] columnIndexes)
    {
        var x = new GeneralMatrix(endRowIndex - startRowIndex + 1, columnIndexes.Length);
        double[][] b = x.Array;
        try
        {
            for (int i = startRowIndex; i <= endRowIndex; i++)
            {
                for (int j = 0; j < columnIndexes.Length; j++)
                {
                    b[i - startRowIndex][j] = _a[i][columnIndexes[j]];
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
        return x;
    }

    /// <summary>
    ///   Get a submatrix.
    /// </summary>
    /// <param name = "r">   Array of row indices.
    /// </param>
    /// <param name = "j0">  Initial column index
    /// </param>
    /// <param name = "j1">  Final column index
    /// </param>
    /// <returns>     A(r(:),j0:j1)
    /// </returns>
    /// <exception cref = "System.IndexOutOfRangeException">   Submatrix indices
    /// </exception>
    public virtual GeneralMatrix GetMatrix(int[] r, int j0, int j1)
    {
        var x = new GeneralMatrix(r.Length, j1 - j0 + 1);
        double[][] b = x.Array;
        try
        {
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = j0; j <= j1; j++)
                {
                    b[i][j - j0] = _a[r[i]][j];
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
        return x;
    }

    /// <summary>
    ///   Set a single element.
    /// </summary>
    /// <param name = "i">   Row index.
    /// </param>
    /// <param name = "j">   Column index.
    /// </param>
    /// <param name = "s">   A(i,j).
    /// </param>
    /// <exception cref = "System.IndexOutOfRangeException">
    /// </exception>
    public virtual void SetElement(int i, int j, double s)
    {
        _a[i][j] = s;
    }

    /// <summary>
    ///   Set a submatrix.
    /// </summary>
    /// <param name = "i0">  Initial row index
    /// </param>
    /// <param name = "i1">  Final row index
    /// </param>
    /// <param name = "j0">  Initial column index
    /// </param>
    /// <param name = "j1">  Final column index
    /// </param>
    /// <param name = "x">   A(startRowIndex:endRowIndex,j0:j1)
    /// </param>
    /// <exception cref = "System.IndexOutOfRangeException">  Submatrix indices
    /// </exception>
    public virtual void SetMatrix(int i0, int i1, int j0, int j1, GeneralMatrix x)
    {
        try
        {
            for (int i = i0; i <= i1; i++)
            {
                for (int j = j0; j <= j1; j++)
                {
                    _a[i][j] = x.GetElement(i - i0, j - j0);
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
    }

    /// <summary>
    ///   Set a submatrix.
    /// </summary>
    /// <param name = "r">   Array of row indices.
    /// </param>
    /// <param name = "c">   Array of column indices.
    /// </param>
    /// <param name = "x">   A(r(:),columnIndexes(:))
    /// </param>
    /// <exception cref = "System.IndexOutOfRangeException">  Submatrix indices
    /// </exception>
    public virtual void SetMatrix(int[] r, int[] c, GeneralMatrix x)
    {
        try
        {
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = 0; j < c.Length; j++)
                {
                    _a[r[i]][c[j]] = x.GetElement(i, j);
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
    }

    /// <summary>
    ///   Set a submatrix.
    /// </summary>
    /// <param name = "r">   Array of row indices.
    /// </param>
    /// <param name = "j0">  Initial column index
    /// </param>
    /// <param name = "j1">  Final column index
    /// </param>
    /// <param name = "x">   A(r(:),j0:j1)
    /// </param>
    /// <exception cref = "System.IndexOutOfRangeException"> Submatrix indices
    /// </exception>
    public virtual void SetMatrix(int[] r, int j0, int j1, GeneralMatrix x)
    {
        try
        {
            for (int i = 0; i < r.Length; i++)
            {
                for (int j = j0; j <= j1; j++)
                {
                    _a[r[i]][j] = x.GetElement(i, j - j0);
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
    }

    /// <summary>
    ///   Set a submatrix.
    /// </summary>
    /// <param name = "i0">  Initial row index
    /// </param>
    /// <param name = "i1">  Final row index
    /// </param>
    /// <param name = "c">   Array of column indices.
    /// </param>
    /// <param name = "x">   A(startRowIndex:endRowIndex,columnIndexes(:))
    /// </param>
    /// <exception cref = "System.IndexOutOfRangeException">  Submatrix indices
    /// </exception>
    public virtual void SetMatrix(int i0, int i1, int[] c, GeneralMatrix x)
    {
        try
        {
            for (int i = i0; i <= i1; i++)
            {
                for (int j = 0; j < c.Length; j++)
                {
                    _a[i][c[j]] = x.GetElement(i - i0, j);
                }
            }
        }
        catch (IndexOutOfRangeException e)
        {
            throw new IndexOutOfRangeException("Submatrix indices", e);
        }
    }

    /// <summary>
    ///   Matrix transpose.
    /// </summary>
    /// <returns>    A'
    /// </returns>
    public virtual GeneralMatrix Transpose()
    {
        var x = new GeneralMatrix(_n, _m);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[j][i] = _a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   One norm
    /// </summary>
    /// <returns>    maximum column sum.
    /// </returns>
    public virtual double Norm1()
    {
        double f = 0;
        for (int j = 0; j < _n; j++)
        {
            double s = 0;
            for (int i = 0; i < _m; i++)
            {
                s += Math.Abs(_a[i][j]);
            }
            f = Math.Max(f, s);
        }
        return f;
    }

    /// <summary>
    ///   Two norm
    /// </summary>
    /// <returns>    maximum singular value.
    /// </returns>
    public virtual double Norm2()
    {
        return (new SingularValueDecomposition(this).Norm2());
    }

    /// <summary>
    ///   Infinity norm
    /// </summary>
    /// <returns>    maximum row sum.
    /// </returns>
    public virtual double NormInf()
    {
        double f = 0;
        for (int i = 0; i < _m; i++)
        {
            double s = 0;
            for (int j = 0; j < _n; j++)
            {
                s += Math.Abs(_a[i][j]);
            }
            f = Math.Max(f, s);
        }
        return f;
    }

    /// <summary>
    ///   Frobenius norm
    /// </summary>
    /// <returns>    sqrt of sum of squares of all elements.
    /// </returns>
    public virtual double NormF()
    {
        double f = 0;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                f = Maths.Hypot(f, _a[i][j]);
            }
        }
        return f;
    }

    /// <summary>
    ///   Unary minus
    /// </summary>
    /// <returns>    -A
    /// </returns>
    public virtual GeneralMatrix UnaryMinus()
    {
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = -_a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   C = A + B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A + B
    /// </returns>
    public virtual GeneralMatrix Add(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = _a[i][j] + b._a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   A = A + B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A + B
    /// </returns>
    public virtual GeneralMatrix AddEquals(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = _a[i][j] + b._a[i][j];
            }
        }
        return this;
    }

    /// <summary>
    ///   C = A - B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A - B
    /// </returns>
    public virtual GeneralMatrix Subtract(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = _a[i][j] - b._a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   A = A - B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A - B
    /// </returns>
    public virtual GeneralMatrix SubtractEquals(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = _a[i][j] - b._a[i][j];
            }
        }
        return this;
    }

    /// <summary>
    ///   Element-by-element multiplication, C = A.*B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A.*B
    /// </returns>
    public virtual GeneralMatrix ArrayMultiply(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = _a[i][j] * b._a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   Element-by-element multiplication in place, A = A.*B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A.*B
    /// </returns>
    public virtual GeneralMatrix ArrayMultiplyEquals(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = _a[i][j] * b._a[i][j];
            }
        }
        return this;
    }

    /// <summary>
    ///   Element-by-element right division, C = A./B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A./B
    /// </returns>
    public virtual GeneralMatrix ArrayRightDivide(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = _a[i][j] / b._a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   Element-by-element right division in place, A = A./B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A./B
    /// </returns>
    public virtual GeneralMatrix ArrayRightDivideEquals(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = _a[i][j] / b._a[i][j];
            }
        }
        return this;
    }

    /// <summary>
    ///   Element-by-element left division, C = A.\B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A.\B
    /// </returns>
    public virtual GeneralMatrix ArrayLeftDivide(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = b._a[i][j] / _a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   Element-by-element left division in place, A = A.\B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     A.\B
    /// </returns>
    public virtual GeneralMatrix ArrayLeftDivideEquals(GeneralMatrix b)
    {
        CheckMatrixDimensions(b);
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = b._a[i][j] / _a[i][j];
            }
        }
        return this;
    }

    /// <summary>
    ///   Multiply a matrix by a scalar, C = s*A
    /// </summary>
    /// <param name = "s">   scalar
    /// </param>
    /// <returns>     s*A
    /// </returns>
    public virtual GeneralMatrix Multiply(double s)
    {
        var x = new GeneralMatrix(_m, _n);
        double[][] c = x.Array;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                c[i][j] = s * _a[i][j];
            }
        }
        return x;
    }

    /// <summary>
    ///   Multiply a matrix by a scalar in place, A = s*A
    /// </summary>
    /// <param name = "s">   scalar
    /// </param>
    /// <returns>     replace A by s*A
    /// </returns>
    public virtual GeneralMatrix MultiplyEquals(double s)
    {
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                _a[i][j] = s * _a[i][j];
            }
        }
        return this;
    }

    /// <summary>
    ///   Linear algebraic matrix multiplication, A * B
    /// </summary>
    /// <param name = "b">   another matrix
    /// </param>
    /// <returns>     Matrix product, A * B
    /// </returns>
    /// <exception cref = "System.ArgumentException">  Matrix inner dimensions must agree.
    /// </exception>
    public virtual GeneralMatrix Multiply(GeneralMatrix b)
    {
        if (b._m != _n)
        {
            throw new ArgumentException("GeneralMatrix inner dimensions must agree.");
        }
        var x = new GeneralMatrix(_m, b._n);
        double[][] c = x.Array;
        var bcolj = new double[_n];
        for (int j = 0; j < b._n; j++)
        {
            for (int k = 0; k < _n; k++)
            {
                bcolj[k] = b._a[k][j];
            }
            for (int i = 0; i < _m; i++)
            {
                double[] arowi = _a[i];
                double s = 0;
                for (int k = 0; k < _n; k++)
                {
                    s += arowi[k] * bcolj[k];
                }
                c[i][j] = s;
            }
        }
        return x;
    }

    /// <summary>
    ///   LU Decomposition
    /// </summary>
    /// <returns>     LUDecomposition
    /// </returns>
    /// <seealso cref = "LUDecomposition">
    /// </seealso>
    public virtual LUDecomposition Lud()
    {
        return new LUDecomposition(this);
    }

    /// <summary>
    ///   QR Decomposition
    /// </summary>
    /// <returns>     QRDecomposition
    /// </returns>
    /// <seealso cref = "QRDecomposition">
    /// </seealso>
    public virtual QRDecomposition Qrd()
    {
        return new QRDecomposition(this);
    }

    /// <summary>
    ///   Cholesky Decomposition
    /// </summary>
    /// <returns>     CholeskyDecomposition
    /// </returns>
    /// <seealso cref = "CholeskyDecomposition">
    /// </seealso>
    public virtual CholeskyDecomposition Chol()
    {
        return new CholeskyDecomposition(this);
    }

    /// <summary>
    ///   Singular Value Decomposition
    /// </summary>
    /// <returns>     SingularValueDecomposition
    /// </returns>
    /// <seealso cref = "SingularValueDecomposition">
    /// </seealso>
    public virtual SingularValueDecomposition Svd()
    {
        return new SingularValueDecomposition(this);
    }

    /// <summary>
    ///   Eigenvalue Decomposition
    /// </summary>
    /// <returns>     EigenvalueDecomposition
    /// </returns>
    /// <seealso cref = "EigenvalueDecomposition">
    /// </seealso>
    public virtual EigenvalueDecomposition Eigen()
    {
        return new EigenvalueDecomposition(this);
    }

    /// <summary>
    ///   Solve A*X = B
    /// </summary>
    /// <param name = "b">   right hand side
    /// </param>
    /// <returns>     solution if A is square, least squares solution otherwise
    /// </returns>
    public virtual GeneralMatrix Solve(GeneralMatrix b)
    {
        return (_m == _n ? (new LUDecomposition(this)).Solve(b) : (new QRDecomposition(this)).Solve(b));
    }

    /// <summary>
    ///   Solve X*A = B, which is also A'*X' = B'
    /// </summary>
    /// <param name = "b">   right hand side
    /// </param>
    /// <returns>     solution if A is square, least squares solution otherwise.
    /// </returns>
    public virtual GeneralMatrix SolveTranspose(GeneralMatrix b)
    {
        return Transpose().Solve(b.Transpose()).Transpose();
    }

    /// <summary>
    ///   Matrix inverse or pseudoinverse
    /// </summary>
    /// <returns>     inverse(A) if A is square, pseudoinverse otherwise.
    /// </returns>
    public virtual GeneralMatrix Inverse()
    {
        return Solve(Identity(_m, _m));
    }

    /// <summary>
    ///   GeneralMatrix determinant
    /// </summary>
    /// <returns>     determinant
    /// </returns>
    public virtual double Determinant()
    {
        return new LUDecomposition(this).Determinant();
    }

    /// <summary>
    ///   GeneralMatrix rank
    /// </summary>
    /// <returns>     effective numerical rank, obtained from SVD.
    /// </returns>
    public virtual int Rank()
    {
        return new SingularValueDecomposition(this).Rank();
    }

    /// <summary>
    ///   Matrix condition (2 norm)
    /// </summary>
    /// <returns>     ratio of largest to smallest singular value.
    /// </returns>
    public virtual double Condition()
    {
        return new SingularValueDecomposition(this).Condition();
    }

    /// <summary>
    ///   Matrix trace.
    /// </summary>
    /// <returns>     sum of the diagonal elements.
    /// </returns>
    public virtual double Trace()
    {
        double t = 0;
        for (int i = 0; i < Math.Min(_m, _n); i++)
        {
            t += _a[i][i];
        }
        return t;
    }

    /// <summary>
    ///   Generate matrix with random elements
    /// </summary>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    /// <returns>     An m-by-n matrix with uniformly distributed random elements.
    /// </returns>
    public static GeneralMatrix Random(int m, int n)
    {
        var random = System.Random.Shared;

        var a = new GeneralMatrix(m, n);
        double[][] x = a.Array;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                x[i][j] = random.NextDouble();
            }
        }
        return a;
    }

    /// <summary>
    ///   Generate matrix with random elements
    /// </summary>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    /// <param name="minValue">lowest value for random number</param>
    /// <param name="maxValue">highest value for random number</param>
    /// <returns>     An m-by-n matrix with uniformly distributed random elements.
    /// </returns>
    public static GeneralMatrix Random(int m, int n, double minValue, double maxValue)
    {
        var random = System.Random.Shared;
        var range = maxValue - minValue;
        var a = new GeneralMatrix(m, n);
        double[][] x = a.Array;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                x[i][j] = random.NextDouble() * range + minValue;
            }
        }
        return a;
    }
    /// <summary>
    ///   Generate matrix with random elements
    /// </summary>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    /// <param name="minValue">lowest value for random number</param>
    /// <param name="maxValue">highest value for random number</param>
    /// <returns>     An m-by-n matrix with uniformly distributed random elements.
    /// </returns>
    public static GeneralMatrix Random(int m, int n, int minValue, int maxValue)
    {
        var random = System.Random.Shared;
        var a = new GeneralMatrix(m, n);
        double[][] x = a.Array;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                x[i][j] = random.Next(minValue, maxValue);
            }
        }
        return a;
    }
    /// <summary>
    ///   Generate identity matrix
    /// </summary>
    /// <param name = "m">   Number of rows.
    /// </param>
    /// <param name = "n">   Number of colums.
    /// </param>
    /// <returns>     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
    /// </returns>
    public static GeneralMatrix Identity(int m, int n)
    {
        var a = new GeneralMatrix(m, n);
        double[][] x = a.Array;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                x[i][j] = (i == j ? 1.0 : 0.0);
            }
        }
        return a;
    }

    #region Operator Overloading

    /// <summary>
    ///   Addition of matrices
    /// </summary>
    /// <param name = "m1"></param>
    /// <param name = "m2"></param>
    /// <returns></returns>
    public static GeneralMatrix operator +(GeneralMatrix m1, GeneralMatrix m2)
    {
        return m1.Add(m2);
    }

    /// <summary>
    ///   Subtraction of matrices
    /// </summary>
    /// <param name = "m1"></param>
    /// <param name = "m2"></param>
    /// <returns></returns>
    public static GeneralMatrix operator -(GeneralMatrix m1, GeneralMatrix m2)
    {
        return m1.Subtract(m2);
    }

    /// <summary>
    ///   Multiplication of matrices
    /// </summary>
    /// <param name = "m1"></param>
    /// <param name = "m2"></param>
    /// <returns></returns>
    public static GeneralMatrix operator *(GeneralMatrix m1, GeneralMatrix m2)
    {
        return m1.Multiply(m2);
    }

    #endregion   //Operator Overloading

    #endregion //  Public Methods

    #region	 Private Methods

    /// <summary>
    ///   Check if size(A) == size(B) *
    /// </summary>
    private void CheckMatrixDimensions(GeneralMatrix b)
    {
        if (b._m != _m || b._n != _n)
        {
            throw new ArgumentException("GeneralMatrix dimensions must agree.");
        }
    }

    #endregion //  Private Methods

    #region Implement IDisposable

    /// <summary>
    ///   Completes disposal. Matrix storage is managed and remains usable for compatibility.
    /// </summary>
    public void Dispose() => GC.SuppressFinalize(this);

    #endregion //  Implement IDisposable

    #region ICloneable Members

    /// <summary>
    ///   Clone the GeneralMatrix object.
    /// </summary>
    public Object Clone()
    {
        return Copy();
    }

    #endregion

    #region ISerializable Members

    /// <summary>
    ///   A method called when serializing this class
    /// </summary>
    /// <param name = "info"></param>
    /// <param name = "context"></param>
    void ISerializable.GetObjectData(SerializationInfo info, StreamingContext context)
    {
    }

    #endregion

    /// <summary>Compares runtime types, dimensions, and each element for exact value equality.</summary>
    /// <remarks>Two null references are equal. NaN values compare equal, as do signed zeros.</remarks>
    /// <param name="m1">The first matrix, or null.</param>
    /// <param name="m2">The second matrix, or null.</param>
    /// <returns>Whether the matrices have the same type, shape, and values.</returns>
    public static bool operator ==(GeneralMatrix? m1, GeneralMatrix? m2)
    {
        return object.Equals(m1, m2);
    }
    /// <summary>Determines whether two matrices differ in runtime type, dimensions, or exact element values.</summary>
    /// <param name="m1">The first matrix, or null.</param>
    /// <param name="m2">The second matrix, or null.</param>
    /// <returns>The negation of the equality operator.</returns>
    public static bool operator !=(GeneralMatrix? m1, GeneralMatrix? m2)
    {
        return !object.Equals(m1, m2);
    }

    /// <summary>Compares this matrix with another matrix using exact element equality.</summary>
    /// <remarks>Runtime types and dimensions must match. NaN values compare equal, as do signed zeros.</remarks>
    /// <param name="other">The matrix to compare, or null.</param>
    /// <returns>Whether the other matrix has the same type, shape, and values.</returns>
    public bool Equals(GeneralMatrix? other)
    {
        if (ReferenceEquals(null, other)) return false;
        if (ReferenceEquals(this, other)) return true;
        if (GetType() != other.GetType()) return false;
        if (_m != other._m || _n != other._n) return false;
        bool result = true;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                if (!_a[i][j].Equals(other._a[i][j]))
                    result = false;
            }
        }
        return result;
    }

    /// <summary>Determines whether an object is an exactly equal matrix of the same runtime type.</summary>
    /// <param name="obj">The object to compare, or null.</param>
    /// <returns>Whether the object is a matrix with the same type, shape, and values.</returns>
    public override bool Equals(object? obj)
    {
        if (ReferenceEquals(null, obj)) return false;
        if (ReferenceEquals(this, obj)) return true;
        return obj is GeneralMatrix matrix && Equals(matrix);
    }

    /// <summary>Computes a hash from the runtime type, dimensions, and element values.</summary>
    /// <remarks>
    /// Changing matrix values can change the hash. Do not mutate a matrix while it is used as a hash key.
    /// The hash is not a persistent identifier and may differ between processes.
    /// </remarks>
    /// <returns>A hash consistent with exact matrix equality.</returns>
    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(GetType());
        hash.Add(_m);
        hash.Add(_n);
        for (int i = 0; i < _m; i++)
            for (int j = 0; j < _n; j++)
                hash.Add(_a[i][j]);
        return hash.ToHashCode();
    }
}
