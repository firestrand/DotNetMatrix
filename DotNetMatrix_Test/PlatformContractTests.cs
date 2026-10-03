using System;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text.Json;
using DotNetMatrix;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace DotNetMatrix_Test;

[TestClass]
public sealed class PlatformContractTests
{
    public TestContext TestContext { get; set; } = null!;

    [TestMethod]
    public void TestProcessRunsOnTheDeclaredNativeArchitecture()
    {
        Architecture observed = RuntimeInformation.ProcessArchitecture;
        TestContext.WriteLine($"Observed OS: {RuntimeInformation.OSDescription}");
        TestContext.WriteLine($"Observed runtime: {RuntimeInformation.FrameworkDescription}");
        TestContext.WriteLine($"Observed process architecture: {observed}");

        string? declared = Environment.GetEnvironmentVariable("EXPECTED_TEST_ARCH");
        bool isCi = string.Equals(Environment.GetEnvironmentVariable("CI"), "true", StringComparison.OrdinalIgnoreCase) ||
            string.Equals(Environment.GetEnvironmentVariable("GITHUB_ACTIONS"), "true", StringComparison.OrdinalIgnoreCase);
        if (isCi)
            Assert.IsFalse(string.IsNullOrWhiteSpace(declared), "CI must explicitly declare EXPECTED_TEST_ARCH for its native test lane.");

        // Local runs still log the actual process; qualification requires a declared CI expectation.
        declared ??= observed.ToString();
        Assert.IsTrue(Enum.TryParse(declared, ignoreCase: true, out Architecture expected) &&
            string.Equals(Enum.GetName(expected), declared, StringComparison.OrdinalIgnoreCase),
            $"EXPECTED_TEST_ARCH must name a supported Architecture value; received '{declared}'.");
        Assert.AreEqual(expected, observed, "Cross-publishing or an installed SDK does not qualify the test process architecture.");
    }

    [TestMethod]
    [DataRow("en-US")]
    [DataRow("fr-FR")]
    [DataRow("tr-TR")]
    [DataRow("ar-SA")]
    public void MatrixWireSchemaIsIndependentOfAmbientCulture(string cultureName)
    {
        CultureInfo originalCulture = CultureInfo.CurrentCulture;
        CultureInfo originalUiCulture = CultureInfo.CurrentUICulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(cultureName);
            CultureInfo.CurrentUICulture = CultureInfo.CurrentCulture;
            var options = new JsonSerializerOptions();
            options.Converters.Add(new MatrixJsonConverter());
            var matrix = new GeneralMatrix(new[] { new[] { 1.5, -2.25, 1e-12 } });
            const string expected = "{\"formatVersion\":1,\"rows\":1,\"columns\":3,\"storageOrder\":\"row-major\",\"values\":[1.5,-2.25,1E-12]}";

            string json = JsonSerializer.Serialize(matrix, options);
            Assert.AreEqual(expected, json, "The persisted schema uses JSON numeric syntax rather than localized number formatting.");
            var restored = JsonSerializer.Deserialize<GeneralMatrix>(expected, options)!;
            CollectionAssert.AreEqual(matrix.RowPackedCopy, restored.RowPackedCopy);
            Assert.ThrowsExactly<JsonException>(() => JsonSerializer.Deserialize<GeneralMatrix>(
                "{\"formatVersion\":1,\"rows\":1,\"columns\":1,\"storageOrder\":\"ROW-MAJOR\",\"values\":[1]}", options));
        }
        finally
        {
            CultureInfo.CurrentCulture = originalCulture;
            CultureInfo.CurrentUICulture = originalUiCulture;
        }
    }
}
