using Microsoft.VisualStudio.TestTools.UnitTesting;

[assembly: DoNotParallelize]

namespace DotNetMatrix_Test;

[TestClass]
public sealed class LegacyHarnessTests
{
    [TestMethod]
    public void LegacyNumericalHarnessReportsNoErrors()
    {
        Assert.AreEqual(0, DotNetMatrix.test.TestMatrix.Run([]));
    }
}
