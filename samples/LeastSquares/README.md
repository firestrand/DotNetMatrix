# Rank-aware least squares

Run from the repository root:

```sh
dotnet restore samples/LeastSquares/LeastSquares.csproj --locked-mode
dotnet run --project samples/LeastSquares/LeastSquares.csproj -c Release --no-restore
```

The analytic examples verify a known intercept/slope fit and duplicated columns.
A deficient design has no unique coefficient vector: the SVD solver returns the
minimum-norm vector at its retained numerical rank. It does not imply that fitted
coefficients are identifiable or statistically meaningful. Residuals describe the
original system; changing the cutoff can discard information and increase them.
No applied helper is added to the library until an actual consumer needs one.
