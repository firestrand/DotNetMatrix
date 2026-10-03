using System;
using System.Numerics;
using BenchmarkDotNet.Attributes;

namespace DotNetMatrix.Benchmarks;

/// <summary>Isolated experiments; these methods do not change the production API.</summary>
[MemoryDiagnoser]
public class MultiplicationCandidates
{
    [Params(32, 128)] public int Size { get; set; }
    private GeneralMatrix left = null!, right = null!, destination = null!;
    private double[] column = null!;

    [GlobalSetup]
    public void Setup()
    {
        left = Fixtures.Create(Size, Size);
        right = Fixtures.Create(Size, Size);
        destination = new GeneralMatrix(Size, Size);
        column = new double[Size];
        AssertCandidates(left, right);
    }

    [Benchmark(Baseline = true)] public GeneralMatrix ExistingMultiply() => left.Multiply(right);
    [Benchmark] public GeneralMatrix CacheBlocked() => Blocked(left, right);
    [Benchmark] public GeneralMatrix ReusedDestinationAndScratch() => Reuse(left, right, destination, column);
    [Benchmark] public GeneralMatrix SimdWithPacking() => Simd(left, right);

    private static GeneralMatrix Reuse(GeneralMatrix a, GeneralMatrix b, GeneralMatrix output, double[] scratch)
    {
        var aa = a.Array; var bb = b.Array; var cc = output.Array;
        for (int j = 0; j < b.ColumnDimension; j++)
        {
            for (int k = 0; k < a.ColumnDimension; k++) scratch[k] = bb[k][j];
            for (int i = 0; i < a.RowDimension; i++)
            {
                double sum = 0;
                for (int k = 0; k < a.ColumnDimension; k++) sum += aa[i][k] * scratch[k];
                cc[i][j] = sum;
            }
        }
        return output;
    }

    private static GeneralMatrix Blocked(GeneralMatrix a, GeneralMatrix b)
    {
        var output = new GeneralMatrix(a.RowDimension, b.ColumnDimension);
        var aa = a.Array; var bb = b.Array; var cc = output.Array;
        const int block = 32;
        for (int ii = 0; ii < a.RowDimension; ii += block)
            for (int jj = 0; jj < b.ColumnDimension; jj += block)
                for (int kk = 0; kk < a.ColumnDimension; kk += block)
                    for (int i = ii; i < Math.Min(ii + block, a.RowDimension); i++)
                        for (int j = jj; j < Math.Min(jj + block, b.ColumnDimension); j++)
                        {
                            double sum = cc[i][j];
                            for (int k = kk; k < Math.Min(kk + block, a.ColumnDimension); k++)
                                sum += aa[i][k] * bb[k][j];
                            cc[i][j] = sum;
                        }
        return output;
    }

    private static GeneralMatrix Simd(GeneralMatrix a, GeneralMatrix b)
    {
        // Include packing and its allocation: this candidate does not assume a cached transpose.
        var packed = b.Transpose().Array;
        var output = new GeneralMatrix(a.RowDimension, b.ColumnDimension);
        var aa = a.Array; var cc = output.Array;
        for (int i = 0; i < a.RowDimension; i++)
            for (int j = 0; j < b.ColumnDimension; j++)
            {
                var sumVector = Vector<double>.Zero;
                int k = 0;
                for (; k <= a.ColumnDimension - Vector<double>.Count; k += Vector<double>.Count)
                    sumVector += new Vector<double>(aa[i], k) * new Vector<double>(packed[j], k);
                double sum = Vector.Sum(sumVector);
                for (; k < a.ColumnDimension; k++) sum += aa[i][k] * packed[j][k];
                cc[i][j] = sum;
            }
        return output;
    }

    public static void Validate()
    {
        foreach (int size in new[] { 1, 3, 16, 32, 64, 128, 192 })
            foreach (var shape in new[] { MatrixShape.Square, MatrixShape.Tall, MatrixShape.Wide })
                AssertCandidates(Fixtures.Create(size, shape), Fixtures.Create(shape == MatrixShape.Wide ? 2 * size : size, size + 1));
        // Cancellation and vastly different magnitudes expose SIMD reassociation.
        AssertCandidates(new GeneralMatrix(new[] { new[] { 1e100, 1, -1e100, 2d } }),
            new GeneralMatrix(new[] { new[] { 1d }, new[] { 1d }, new[] { 1d }, new[] { 1d } }), cancellation: true);
    }

    private static void AssertCandidates(GeneralMatrix a, GeneralMatrix b, bool cancellation = false)
    {
        var aCopy = a.Copy(); var bCopy = b.Copy();
        var expected = a.Multiply(b);
        var output = new GeneralMatrix(a.RowDimension, b.ColumnDimension);
        var scratch = new double[a.ColumnDimension];
        Check(expected, Blocked(a, b), a, b, cancellation);
        Check(expected, Reuse(a, b, output, scratch), a, b, cancellation);
        output.Array[0][0] = double.NaN;
        Check(expected, Reuse(a, b, output, scratch), a, b, cancellation);
        Check(expected, Simd(a, b), a, b, cancellation);
        if (!a.Equals(aCopy) || !b.Equals(bCopy)) throw new InvalidOperationException("Prototype changed input storage.");
    }

    private static void Check(GeneralMatrix expected, GeneralMatrix actual, GeneralMatrix a, GeneralMatrix b, bool cancellation)
    {
        double error = expected.Subtract(actual).NormF();
        double bound = 1e-12 * a.NormF() * b.NormF();
        if (!double.IsFinite(error) || error > bound)
            throw new InvalidOperationException($"Candidate numerical error {error:R} exceeds {bound:R}.");
        if (cancellation && error != 0)
            Console.WriteLine($"Cancellation fixture: reassociated result differs by {error:R}; normwise bound {bound:R}.");
    }
}
