using System;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace DotNetMatrix;

/// <summary>Opt-in JSON persistence for finite matrices using the version 1 row-major schema.</summary>
/// <remarks>
/// Register this converter in <see cref="JsonSerializerOptions.Converters"/>. The schema contains
/// formatVersion, rows, columns, storageOrder, and values. Unknown, duplicate, or missing fields
/// and null matrices are rejected. Existing legacy serialization and matrix aliases are unchanged.
/// </remarks>
public sealed class MatrixJsonConverter : JsonConverter<GeneralMatrix>
{
    private readonly int _maxElements;

    /// <summary>Creates a converter with a bound on numeric elements and allocated jagged rows.</summary>
    /// <param name="maxElements">Positive allocation limit; the default is one million.</param>
    /// <exception cref="ArgumentOutOfRangeException">The limit is not positive.</exception>
    /// <remarks>
    /// The row count also cannot exceed this limit, including matrices with zero columns.
    /// A zero-row matrix may have any nonnegative Int32 column count since it allocates no rows.
    /// Caller-controlled input byte limits remain the responsibility of the transport or serializer.
    /// </remarks>
    public MatrixJsonConverter(int maxElements = 1_000_000)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxElements);
        _maxElements = maxElements;
    }

    /// <summary>Handles null explicitly so the finite-matrix schema cannot be bypassed.</summary>
    public override bool HandleNull => true;

    /// <summary>Reads and validates a complete version 1 matrix without allocating from unchecked dimensions.</summary>
    /// <exception cref="JsonException">The payload does not meet the schema or allocation limits.</exception>
    public override GeneralMatrix Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
    {
        if (reader.TokenType != JsonTokenType.StartObject)
        {
            throw new JsonException("Expected a matrix object.");
        }

        int fields = 0, rows = 0, columns = 0, valueCount = 0;
        Utf8JsonReader valuesReader = default;
        while (reader.Read())
        {
            if (reader.TokenType == JsonTokenType.EndObject)
            {
                if (fields != 31)
                {
                    throw new JsonException("All five matrix fields are required.");
                }

                int elementCount = ValidateShape(rows, columns);
                if (valueCount != elementCount)
                {
                    throw new JsonException("Values count must match rows times columns.");
                }

                var matrix = new GeneralMatrix(rows, columns);
                for (int row = 0; row < rows; row++)
                {
                    for (int column = 0; column < columns; column++)
                    {
                        valuesReader.Read();
                        matrix.SetElement(row, column, valuesReader.GetDouble());
                    }
                }
                return matrix;
            }

            if (reader.TokenType != JsonTokenType.PropertyName)
            {
                throw new JsonException("Expected a matrix field.");
            }

            int field = reader.ValueTextEquals("formatVersion") ? 1
                : reader.ValueTextEquals("rows") ? 2
                : reader.ValueTextEquals("columns") ? 4
                : reader.ValueTextEquals("storageOrder") ? 8
                : reader.ValueTextEquals("values") ? 16
                : throw new JsonException("Unknown matrix field.");
            if ((fields & field) != 0)
            {
                throw new JsonException("Duplicate matrix field.");
            }
            fields |= field;
            if (!reader.Read())
            {
                throw new JsonException("Missing matrix field value.");
            }

            switch (field)
            {
                case 1:
                    if (ReadInteger(ref reader) != 1)
                    {
                        throw new JsonException("Unsupported matrix format version.");
                    }
                    break;
                case 2:
                    rows = ReadInteger(ref reader);
                    if (rows < 0 || rows > _maxElements)
                    {
                        throw new JsonException("Matrix row count exceeds the allocation limit.");
                    }
                    break;
                case 4:
                    columns = ReadInteger(ref reader);
                    if (columns < 0)
                    {
                        throw new JsonException("Matrix column count must be nonnegative.");
                    }
                    break;
                case 8:
                    if (reader.TokenType != JsonTokenType.String || !reader.ValueTextEquals("row-major"))
                    {
                        throw new JsonException("Unsupported matrix storage order.");
                    }
                    break;
                case 16:
                    if (reader.TokenType != JsonTokenType.StartArray)
                    {
                        throw new JsonException("Matrix values must be an array.");
                    }
                    valuesReader = reader;
                    while (reader.Read() && reader.TokenType != JsonTokenType.EndArray)
                    {
                        if (valueCount == _maxElements)
                        {
                            throw new JsonException("Matrix values exceed the element limit.");
                        }
                        if (reader.TokenType != JsonTokenType.Number || !reader.TryGetDouble(out double value) || !double.IsFinite(value))
                        {
                            throw new JsonException("Matrix values must be finite JSON numbers.");
                        }
                        valueCount++;
                    }
                    if (reader.TokenType != JsonTokenType.EndArray)
                    {
                        throw new JsonException("Incomplete matrix values array.");
                    }
                    break;
            }

            if ((fields & 6) == 6)
            {
                ValidateShape(rows, columns);
            }
        }

        throw new JsonException("Incomplete matrix object.");
    }

    /// <summary>Writes a finite matrix in row-major order without changing its storage.</summary>
    /// <exception cref="JsonException">The matrix is null, contains nonfinite values, or exceeds allocation limits.</exception>
    /// <remarks>Callers must prevent concurrent modification of the matrix and its borrowed arrays.</remarks>
    public override void Write(Utf8JsonWriter writer, GeneralMatrix value, JsonSerializerOptions options)
    {
        if (value is null)
        {
            throw new JsonException("Null is not a matrix payload.");
        }

        int rows = value.RowDimension, columns = value.ColumnDimension;
        ValidateShape(rows, columns);
        for (int row = 0; row < rows; row++)
        {
            for (int column = 0; column < columns; column++)
            {
                if (!double.IsFinite(value.GetElement(row, column)))
                {
                    throw new JsonException("Matrix values must be finite JSON numbers.");
                }
            }
        }

        writer.WriteStartObject();
        writer.WriteNumber("formatVersion", 1);
        writer.WriteNumber("rows", rows);
        writer.WriteNumber("columns", columns);
        writer.WriteString("storageOrder", "row-major");
        writer.WriteStartArray("values");
        for (int row = 0; row < rows; row++)
        {
            for (int column = 0; column < columns; column++)
            {
                writer.WriteNumberValue(value.GetElement(row, column));
            }
        }
        writer.WriteEndArray();
        writer.WriteEndObject();
    }

    private int ValidateShape(int rows, int columns)
    {
        long elements = checked((long)rows * columns);
        if (rows < 0 || columns < 0 || rows > _maxElements || elements > _maxElements)
        {
            throw new JsonException("Matrix dimensions exceed the allocation limit.");
        }
        return (int)elements;
    }

    private static int ReadInteger(ref Utf8JsonReader reader)
    {
        if (reader.TokenType != JsonTokenType.Number || !reader.TryGetInt32(out int value))
        {
            throw new JsonException("Matrix metadata must contain Int32 numbers.");
        }
        return value;
    }
}
