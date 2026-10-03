using System;

namespace DotNetMatrix;

internal static class NumericalGuard
{
    internal static void Matrix(GeneralMatrix matrix, string parameterName, bool square = false, bool tall = false, bool nonempty = true)
    {
        ArgumentNullException.ThrowIfNull(matrix, parameterName);
        int rows = matrix.RowDimension, columns = matrix.ColumnDimension;
        if (rows < 0 || columns < 0 || (nonempty && (rows == 0 || columns == 0)) ||
            (square && rows != columns) || (tall && rows < columns))
        {
            throw new ArgumentException("Matrix dimensions are not supported by this operation.", parameterName);
        }
        var values = matrix.Array;
        if (values.Length != rows)
        {
            throw new ArgumentException("Matrix storage does not match its dimensions.", parameterName);
        }
        for (int i = 0; i < rows; i++)
        {
            if (values[i] is null || values[i].Length != columns)
            {
                throw new ArgumentException("Matrix storage does not match its dimensions.", parameterName);
            }
            for (int j = 0; j < columns; j++)
            {
                if (!double.IsFinite(values[i][j]))
                {
                    throw new ArgumentException("Numerical operations require finite matrix elements.", parameterName);
                }
            }
        }
    }

    internal static int Iterations(int maxIterations)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxIterations);
        return maxIterations;
    }
}
