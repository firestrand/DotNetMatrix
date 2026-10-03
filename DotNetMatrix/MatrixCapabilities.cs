using System;
using System.Diagnostics.CodeAnalysis;

namespace DotNetMatrix;

public partial class GeneralMatrix
{
    /// <summary>Gets or sets an element in the existing mutable storage.</summary>
    public double this[int row, int column]
    {
        get => GetElement(row, column);
        set => SetElement(row, column, value);
    }

    /// <summary>Returns an independent 1-by-columns copy of a row.</summary>
    public GeneralMatrix GetRow(int row)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(row);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(row, _m);
        var result = new GeneralMatrix(1, _n);
        _a[row].AsSpan(0, _n).CopyTo(result._a[0]);
        return result;
    }

    /// <summary>Returns an independent rows-by-1 copy of a column.</summary>
    public GeneralMatrix GetColumn(int column)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(column);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(column, _n);
        var result = new GeneralMatrix(_m, 1);
        for (int i = 0; i < _m; i++) result._a[i][0] = _a[i][column];
        return result;
    }

    /// <summary>Creates a square matrix with copied diagonal values.</summary>
    public static GeneralMatrix Diagonal(params double[] values)
    {
        ArgumentNullException.ThrowIfNull(values);
        var result = new GeneralMatrix(values.Length, values.Length);
        for (int i = 0; i < values.Length; i++) result._a[i][i] = values[i];
        return result;
    }

    /// <summary>Stacks rows in independently owned storage; column counts must agree.</summary>
    public GeneralMatrix ConcatenateRows(GeneralMatrix other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (_n != other._n) throw new ArgumentException("Column counts must agree.", nameof(other));
        var result = new GeneralMatrix(checked(_m + other._m), _n);
        for (int i = 0; i < _m; i++) _a[i].AsSpan(0, _n).CopyTo(result._a[i]);
        for (int i = 0; i < other._m; i++) other._a[i].AsSpan(0, _n).CopyTo(result._a[_m + i]);
        return result;
    }

    /// <summary>Joins columns in independently owned storage; row counts must agree.</summary>
    public GeneralMatrix ConcatenateColumns(GeneralMatrix other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (_m != other._m) throw new ArgumentException("Row counts must agree.", nameof(other));
        var result = new GeneralMatrix(_m, checked(_n + other._n));
        for (int i = 0; i < _m; i++)
        {
            _a[i].AsSpan(0, _n).CopyTo(result._a[i]);
            other._a[i].AsSpan(0, other._n).CopyTo(result._a[i].AsSpan(_n));
        }
        return result;
    }

    /// <summary>Copies to caller storage, row-major by default; trailing buffer elements are untouched.</summary>
    /// <remarks>The destination must not overlap borrowed matrix storage.</remarks>
    public void CopyTo(Span<double> destination, bool columnMajor = false)
    {
        if ((long)_m * _n > destination.Length) throw new ArgumentException("Destination is too short.", nameof(destination));
        int index = 0;
        if (columnMajor)
        {
            for (int j = 0; j < _n; j++)
                for (int i = 0; i < _m; i++) destination[index++] = _a[i][j];
        }
        else
        {
            for (int i = 0; i < _m; i++)
                for (int j = 0; j < _n; j++) destination[index++] = _a[i][j];
        }
    }

    /// <summary>Tests |a-b| ≤ absoluteTolerance + relativeTolerance*max(|a|,|b|) elementwise.</summary>
    /// <remarks>NaNs never compare approximately equal. Equal infinities and signed zeros do.
    /// Finite comparisons are scaled to avoid intermediate overflow. This is not hash equality.</remarks>
    public bool ApproximatelyEquals(GeneralMatrix? other, double absoluteTolerance = 1e-12, double relativeTolerance = 1e-12)
    {
        ValidateTolerance(absoluteTolerance, nameof(absoluteTolerance));
        ValidateTolerance(relativeTolerance, nameof(relativeTolerance));
        if (other is null || _m != other._m || _n != other._n) return false;
        for (int i = 0; i < _m; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                double a = _a[i][j], b = other._a[i][j];
                if (double.IsNaN(a) || double.IsNaN(b)) return false;
                if (a == b) continue;
                if (!double.IsFinite(a) || !double.IsFinite(b)) return false;
                double scale = Math.Max(Math.Abs(a), Math.Abs(b));
                // Direct difference preserves close-value precision when it cannot overflow.
                double difference = Math.Abs(a - b);
                double scaledDifference = double.IsFinite(difference) ? difference / scale : Math.Abs(a / scale - b / scale);
                if (scaledDifference > absoluteTolerance / scale + relativeTolerance) return false;
            }
        }
        return true;
    }

    /// <summary>Uses a caller-owned random generator. Streams are runtime-dependent.</summary>
    public static GeneralMatrix Random(int rows, int columns, System.Random random)
    {
        ArgumentNullException.ThrowIfNull(random);
        ArgumentOutOfRangeException.ThrowIfNegative(rows);
        ArgumentOutOfRangeException.ThrowIfNegative(columns);
        var result = new GeneralMatrix(rows, columns);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++) result._a[i][j] = random.NextDouble();
        return result;
    }

    /// <summary>Computes an SVD Moore–Penrose inverse of the retained-rank matrix.</summary>
    /// <param name="relativeTolerance">Finite nonnegative cutoff relative to largest singular value;
    /// values exactly at the cutoff are discarded. Default is max(rows,columns)*2^-52.</param>
    public GeneralMatrix PseudoInverse(double? relativeTolerance = null)
    {
        double tolerance = Cutoff(relativeTolerance);
        var svd = new SingularValueDecomposition(this);
        return PseudoInverse(svd, tolerance);
    }

    /// <summary>Returns the minimum-norm least-squares solution at the retained SVD rank.</summary>
    public GeneralMatrix SolveMinimumNorm(GeneralMatrix rightHandSide, double? relativeTolerance = null)
        => SolveMinimumNormWithDiagnostics(rightHandSide, relativeTolerance).Solution;

    /// <summary>Computes minimum-norm solution and diagnostics using one SVD.</summary>
    /// <remarks>A failed convergence throws; successful results have Converged=true.
    /// Condition number, when requested, is for the original matrix, not its truncation.</remarks>
    public MatrixSolveResult SolveMinimumNormWithDiagnostics(GeneralMatrix rightHandSide, double? relativeTolerance = null, bool includeCondition = false)
    {
        ValidateRightHandSide(rightHandSide);
        double tolerance = Cutoff(relativeTolerance);
        var svd = new SingularValueDecomposition(this);
        var solution = new GeneralMatrix(_n, rightHandSide._n);
        double[][] u = svd.GetU().Array, v = svd.GetV().Array;
        var columnScales = new double[rightHandSide._n];
        for (int j = 0; j < rightHandSide._n; j++)
            for (int i = 0; i < _m; i++) columnScales[j] = Math.Max(columnScales[j], Math.Abs(rightHandSide._a[i][j]));
        double threshold = tolerance * svd.SingularValues[0];
        int rank = 0;
        for (int k = 0; k < svd.SingularValues.Length; k++)
        {
            double sigma = svd.SingularValues[k];
            if (sigma <= threshold) continue;
            rank++;
            for (int j = 0; j < rightHandSide._n; j++)
            {
                double scale = columnScales[j];
                if (scale == 0) continue;
                double projection = 0;
                for (int i = 0; i < _m; i++) projection += u[i][k] * (rightHandSide._a[i][j] / scale);
                // Apply the scale/sigma ratio in exponent space: either intermediate
                // quotient can overflow even when the final coefficient is representable.
                int scaleExponent = Math.ILogB(scale), sigmaExponent = Math.ILogB(sigma);
                double coefficient = Math.ScaleB(projection *
                    (Math.ScaleB(scale, -scaleExponent) / Math.ScaleB(sigma, -sigmaExponent)),
                    scaleExponent - sigmaExponent);
                for (int i = 0; i < _n; i++) solution._a[i][j] += v[i][k] * coefficient;
            }
        }
        return new MatrixSolveResult(solution, Multiply(solution).Subtract(rightHandSide).NormF(), rank,
            includeCondition ? svd.Condition() : null);
    }

    /// <summary>Uses the existing LU/QR solve and computes the Frobenius residual.</summary>
    /// <remarks>Rank is full column rank when the ordinary solve succeeds. Optional conditioning
    /// requires a separate SVD; it is not computed by default.</remarks>
    public MatrixSolveResult SolveWithDiagnostics(GeneralMatrix rightHandSide, bool includeCondition = false)
    {
        ValidateRightHandSide(rightHandSide);
        var solution = Solve(rightHandSide);
        return new MatrixSolveResult(solution, Multiply(solution).Subtract(rightHandSide).NormF(), _n,
            includeCondition ? Condition() : null);
    }

    /// <summary>Returns false for singular/rank-deficient systems. Invalid inputs remain errors.</summary>
    public bool TrySolve(GeneralMatrix rightHandSide, [NotNullWhen(true)] out GeneralMatrix? solution)
    {
        ValidateRightHandSide(rightHandSide);
        solution = null;
        if (_m == _n)
        {
            var lu = new LUDecomposition(this);
            if (!lu.IsNonSingular) return false;
            solution = lu.Solve(rightHandSide);
        }
        else
        {
            var qr = new QRDecomposition(this);
            if (!qr.FullRank) return false;
            solution = qr.Solve(rightHandSide);
        }
        return true;
    }

    private GeneralMatrix PseudoInverse(SingularValueDecomposition svd, double tolerance)
    {
        var result = new GeneralMatrix(_n, _m);
        double[] singularValues = svd.SingularValues;
        double threshold = tolerance * singularValues[0];
        double[][] u = svd.GetU().Array, v = svd.GetV().Array;
        for (int k = 0; k < Math.Min(_m, _n); k++)
        {
            if (singularValues[k] <= threshold) continue;
            for (int i = 0; i < _n; i++)
                for (int j = 0; j < _m; j++) result._a[i][j] += v[i][k] * (u[j][k] / singularValues[k]);
        }
        return result;
    }

    private double Cutoff(double? relativeTolerance)
    {
        double tolerance = relativeTolerance ?? Math.Max(_m, _n) * Math.ScaleB(1.0, -52);
        ValidateTolerance(tolerance, nameof(relativeTolerance));
        return tolerance;
    }

    private static void ValidateTolerance(double value, string name)
    {
        if (!double.IsFinite(value) || value < 0) throw new ArgumentOutOfRangeException(name, "Tolerance must be finite and nonnegative.");
    }

    private void ValidateRightHandSide(GeneralMatrix rightHandSide)
    {
        ArgumentNullException.ThrowIfNull(rightHandSide);
        if (rightHandSide._m != _m) throw new ArgumentException("Right-hand side row count must agree.", nameof(rightHandSide));
        NumericalGuard.Matrix(rightHandSide, nameof(rightHandSide), nonempty: false);
    }
}

/// <summary>Solution and diagnostics. The solution has independently owned mutable storage.</summary>
public sealed class MatrixSolveResult
{
    internal MatrixSolveResult(GeneralMatrix solution, double residualNorm, int numericalRank, double? conditionNumber)
    {
        Solution = solution;
        ResidualNorm = residualNorm;
        NumericalRank = numericalRank;
        ConditionNumber = conditionNumber;
    }

    /// <summary>Computed solution.</summary>
    public GeneralMatrix Solution { get; }
    /// <summary>Frobenius norm of A*X-B for the original system.</summary>
    public double ResidualNorm { get; }
    /// <summary>Numerical rank using the solving algorithm's documented policy.</summary>
    public int NumericalRank { get; }
    /// <summary>Original matrix two-norm condition number, or null when not requested.</summary>
    public double? ConditionNumber { get; }
    /// <summary>True on success; convergence failures throw rather than return partial results.</summary>
    public bool Converged => true;
}
