using System;
using System.IO;
using System.Text.Json;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class MatrixJsonTests
{
    private static JsonSerializerOptions Options(int maxElements = 1_000_000)
    {
        var options = new JsonSerializerOptions();
        options.Converters.Add(new MatrixJsonConverter(maxElements));
        return options;
    }

    private static string Payload(string rows, string columns, string values) =>
        $"{{\"formatVersion\":1,\"rows\":{rows},\"columns\":{columns},\"storageOrder\":\"row-major\",\"values\":{values}}}";

    [TestMethod]
    public void RectangularRoundtripUsesExplicitRowMajorSchemaAndIndependentStorage()
    {
        var matrix = new GeneralMatrix(new[] { new[] { 1.0, -2, 3.5 }, new[] { 4.0, 0, 6 } });
        string json = JsonSerializer.Serialize(matrix, Options());
        using var document = JsonDocument.Parse(json);
        var root = document.RootElement;
        Assert.AreEqual(5, root.GetPropertyCount());
        Assert.AreEqual(1, root.GetProperty("formatVersion").GetInt32());
        Assert.AreEqual(2, root.GetProperty("rows").GetInt32());
        Assert.AreEqual(3, root.GetProperty("columns").GetInt32());
        Assert.AreEqual("row-major", root.GetProperty("storageOrder").GetString());
        Assert.AreEqual("[1,-2,3.5,4,0,6]", root.GetProperty("values").GetRawText());
        var restored = JsonSerializer.Deserialize<GeneralMatrix>(json, Options())!;
        Assert.AreEqual(2, restored.RowDimension);
        Assert.AreEqual(3, restored.ColumnDimension);
        CollectionAssert.AreEqual(matrix.RowPackedCopy, restored.RowPackedCopy);
        restored.SetElement(0, 0, 99);
        Assert.AreEqual(1.0, matrix.GetElement(0, 0));
    }

    [TestMethod]
    [DataRow(0, 0)]
    [DataRow(0, 4)]
    [DataRow(3, 0)]
    public void EmptyMatricesPreserveBothDimensions(int rows, int columns)
    {
        var matrix = new GeneralMatrix(rows, columns);
        var restored = JsonSerializer.Deserialize<GeneralMatrix>(JsonSerializer.Serialize(matrix, Options()), Options())!;
        Assert.AreEqual(rows, restored.RowDimension);
        Assert.AreEqual(columns, restored.ColumnDimension);
        Assert.AreEqual(0, restored.RowPackedCopy.Length);
    }

    [TestMethod]
    public void PropertyOrderAndSerializerNamingPoliciesDoNotChangeTheSchema()
    {
        const string json = "{\"values\":[2,3],\"columns\":2,\"storageOrder\":\"row-major\",\"rows\":1,\"formatVersion\":1}";
        var options = Options();
        options.PropertyNamingPolicy = JsonNamingPolicy.SnakeCaseUpper;
        var matrix = JsonSerializer.Deserialize<GeneralMatrix>(json, options)!;
        CollectionAssert.AreEqual(new[] { 2.0, 3 }, matrix.RowPackedCopy);
        StringAssert.Contains(JsonSerializer.Serialize(matrix, options), "\"formatVersion\"");
    }

    [TestMethod]
    [DataRow("{}")]
    [DataRow("[]")]
    [DataRow("null")]
    [DataRow("42")]
    [DataRow("{\"formatVersion\":2,\"rows\":0,\"columns\":0,\"storageOrder\":\"row-major\",\"values\":[]}")]
    [DataRow("{\"formatVersion\":1,\"rows\":0,\"columns\":0,\"storageOrder\":\"column-major\",\"values\":[]}")]
    [DataRow("{\"formatVersion\":1,\"rows\":0,\"columns\":0,\"storageOrder\":\"row-major\",\"values\":[],\"extra\":true}")]
    public void UnknownOrIncompleteSchemasAreRejected(string json)
    {
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(json, Options()));
    }

    [TestMethod]
    public void EveryRequiredFieldMustAppearExactlyOnce()
    {
        string[] fields = { "\"formatVersion\":1", "\"rows\":0", "\"columns\":0", "\"storageOrder\":\"row-major\"", "\"values\":[]" };
        for (int missing = 0; missing < fields.Length; missing++)
        {
            string omitted = "{" + string.Join(',', System.Array.FindAll(fields, field => field != fields[missing])) + "}";
            Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(omitted, Options()));
            string duplicate = "{" + string.Join(',', fields) + "," + fields[missing] + "}";
            Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(duplicate, Options()));
        }
        const string escapedDuplicate = "{\"formatVersion\":1,\"rows\":0,\"ro\\u0077s\":0,\"columns\":0,\"storageOrder\":\"row-major\",\"values\":[]}";
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(escapedDuplicate, Options()));
    }

    [TestMethod]
    [DataRow("-1", "0")]
    [DataRow("0", "-1")]
    [DataRow("1.5", "0")]
    [DataRow("\"1\"", "0")]
    [DataRow("2147483648", "0")]
    [DataRow("2147483647", "2147483647")]
    [DataRow("2147483647", "0")]
    public void InvalidOrExcessiveDimensionsAreRejectedWithoutMatrixAllocation(string rows, string columns)
    {
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(Payload(rows, columns, "[]"), Options()));
    }

    [TestMethod]
    [DataRow("[]")]
    [DataRow("[1,2]")]
    [DataRow("null")]
    [DataRow("{}")]
    [DataRow("[\"1\"]")]
    [DataRow("[null]")]
    [DataRow("[true]")]
    [DataRow("[[]]")]
    [DataRow("[{}]")]
    [DataRow("[1e999]")]
    [DataRow("[\"NaN\"]")]
    public void InvalidCountsAndNonnumericOrNonfiniteValuesAreRejected(string values)
    {
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(Payload("1", "1", values), Options()));
    }

    [TestMethod]
    public void ElementAndJaggedRowLimitsAreEnforcedOnReadAndWrite()
    {
        var options = Options(2);
        var allowed = JsonSerializer.Deserialize<GeneralMatrix>(Payload("1", "2", "[3,4]"), options)!;
        CollectionAssert.AreEqual(new[] { 3.0, 4 }, allowed.RowPackedCopy);
        Assert.AreEqual(Payload("1", "2", "[3,4]"), JsonSerializer.Serialize(allowed, options));
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(Payload("1", "3", "[1,2,3]"), options));
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>("{\"values\":[1,2,3],\"formatVersion\":1,\"rows\":1,\"columns\":3,\"storageOrder\":\"row-major\"}", options));
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(Payload("3", "0", "[]"), options));
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Serialize(new GeneralMatrix(1, 3), options));
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Serialize(new GeneralMatrix(3, 0), options));
        var zeroRows = JsonSerializer.Deserialize<GeneralMatrix>(Payload("0", "2147483647", "[]"), options)!;
        Assert.AreEqual(int.MaxValue, zeroRows.ColumnDimension);
    }

    [TestMethod]
    [DataRow(0)]
    [DataRow(-1)]
    public void ConfiguredElementLimitMustBePositive(int limit)
    {
        Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new MatrixJsonConverter(limit));
    }

    [TestMethod]
    public void FiniteExtremesAndSignedZeroRoundtrip()
    {
        var matrix = new GeneralMatrix(new[] { new[] { double.MaxValue, double.MinValue, double.Epsilon, -0.0 } });
        var restored = JsonSerializer.Deserialize<GeneralMatrix>(JsonSerializer.Serialize(matrix, Options()), Options())!;
        CollectionAssert.AreEqual(matrix.RowPackedCopy, restored.RowPackedCopy);
        Assert.AreEqual(BitConverter.DoubleToInt64Bits(-0.0), BitConverter.DoubleToInt64Bits(restored.GetElement(0, 3)));
    }

    [TestMethod]
    [DataRow(double.NaN)]
    [DataRow(double.PositiveInfinity)]
    [DataRow(double.NegativeInfinity)]
    public void SerializationRejectsNonfiniteValuesBeforeWritingPayload(double invalid)
    {
        var converter = new MatrixJsonConverter();
        using var stream = new MemoryStream();
        using var writer = new Utf8JsonWriter(stream);
        Assert.ThrowsExactly<JsonException>(() => converter.Write(writer, new GeneralMatrix(1, 1, invalid), Options()));
        writer.Flush();
        Assert.AreEqual(0L, stream.Length);
    }

    [TestMethod]
    public void DirectConverterCallsRejectNullAndIncompleteReaders()
    {
        Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Serialize<GeneralMatrix>(null!, Options()));
        var reader = new Utf8JsonReader("{\"rows\":0}"u8);
        var converter = new MatrixJsonConverter();
        try
        {
            converter.Read(ref reader, typeof(GeneralMatrix), Options());
            Assert.Fail("An unstarted reader should be rejected.");
        }
        catch (JsonException) { }
    }
}
