using System;
using BenchmarkDotNet.Attributes;

namespace DotNetMatrix.Benchmarks;

public enum MatrixShape { Square, Tall, Wide }

internal static class Fixtures
{
    internal static GeneralMatrix Create(int rows, int columns)
    {
        var result = new GeneralMatrix(rows, columns);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                result.SetElement(i, j, Math.Sin((i + 1) * (j + 2)) + (i == j ? 2 : 0));
        return result;
    }

    internal static GeneralMatrix Create(int size, MatrixShape shape) => shape switch
    {
        MatrixShape.Tall => Create(2 * size, size),
        MatrixShape.Wide => Create(size, 2 * size),
        _ => Create(size, size)
    };
}

[MemoryDiagnoser]
public class MatrixOperations
{
    [Params(16, 64, 192)] public int Size { get; set; }
    [Params(MatrixShape.Square, MatrixShape.Tall, MatrixShape.Wide)] public MatrixShape Shape { get; set; }
    private GeneralMatrix input = null!, right = null!;

    [GlobalSetup]
    public void Setup()
    {
        input = Fixtures.Create(Size, Shape);
        right = Fixtures.Create(input.ColumnDimension, Size);
    }

    [Benchmark] public GeneralMatrix Multiply() => input.Multiply(right);
    [Benchmark] public GeneralMatrix Transpose() => input.Transpose();
    [Benchmark] public GeneralMatrix Copy() => input.Copy();
    [Benchmark] public double NormF() => input.NormF();
    [Benchmark] public double Norm1() => input.Norm1();
    [Benchmark] public double NormInfinity() => input.NormInf();
}

[MemoryDiagnoser]
public class TallSquareFactorization
{
    [Params(16, 64, 192)] public int Size { get; set; }
    [Params(MatrixShape.Square, MatrixShape.Tall)] public MatrixShape Shape { get; set; }
    private GeneralMatrix input = null!, right = null!;
    private QRDecomposition qr = null!;

    [GlobalSetup]
    public void Setup()
    {
        input = Fixtures.Create(Size, Shape);
        right = Fixtures.Create(input.RowDimension, 4);
        qr = new QRDecomposition(input);
    }

    [Benchmark] public LUDecomposition FactorLU() => new(input);
    [Benchmark] public QRDecomposition FactorQR() => new(input);
    [Benchmark] public GeneralMatrix RepeatedQrSolve() => qr.Solve(right);
}

[MemoryDiagnoser]
public class SquareFactorization
{
    [Params(16, 64, 192)] public int Size { get; set; }
    private GeneralMatrix input = null!, positiveDefinite = null!, right = null!;
    private LUDecomposition lu = null!;
    private CholeskyDecomposition chol = null!;

    [GlobalSetup]
    public void Setup()
    {
        input = Fixtures.Create(Size, Size);
        positiveDefinite = input.Transpose().Multiply(input).Add(GeneralMatrix.Identity(Size, Size));
        right = Fixtures.Create(Size, 4);
        lu = new LUDecomposition(input);
        chol = new CholeskyDecomposition(positiveDefinite);
    }

    [Benchmark] public CholeskyDecomposition FactorCholesky() => new(positiveDefinite);
    [Benchmark] public EigenvalueDecomposition FactorSymmetricEigen() => new(positiveDefinite);
    [Benchmark] public EigenvalueDecomposition FactorNonsymmetricEigen() => new(input);
    [Benchmark] public GeneralMatrix RepeatedLuSolve() => lu.Solve(right);
    [Benchmark] public GeneralMatrix RepeatedCholeskySolve() => chol.Solve(right);
}

[MemoryDiagnoser]
public class SingularValueFactorization
{
    [Params(16, 64, 192)] public int Size { get; set; }
    [Params(MatrixShape.Square, MatrixShape.Tall, MatrixShape.Wide)] public MatrixShape Shape { get; set; }
    private GeneralMatrix input = null!, right = null!, pseudoinverse = null!;

    [GlobalSetup]
    public void Setup()
    {
        input = Fixtures.Create(Size, Shape);
        right = Fixtures.Create(input.RowDimension, 4);
        pseudoinverse = input.PseudoInverse();
    }

    [Benchmark] public SingularValueDecomposition FactorSVD() => new(input);
    // Isolates repeated application from the factorization and pseudoinverse setup.
    [Benchmark] public GeneralMatrix RepeatedMinimumNormSolve() => pseudoinverse.Multiply(right);
}
