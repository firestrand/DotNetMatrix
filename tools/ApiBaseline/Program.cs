using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata;
using System.Reflection.Metadata.Ecma335;
using System.Reflection.PortableExecutable;
using System.Text.Json;

var commands = new HashSet<string>(StringComparer.Ordinal) { "api", "coverage-types", "source-link" };
if (args.Length != 2 || !commands.Contains(args[0]))
{
    Console.Error.WriteLine("Usage: ApiBaseline api|coverage-types|source-link ASSEMBLY_PATH");
    return 2;
}

if (string.Equals(args[0], "source-link", StringComparison.Ordinal))
{
    using FileStream pdb = File.OpenRead(Path.ChangeExtension(args[1], ".pdb"));
    using MetadataReaderProvider provider = MetadataReaderProvider.FromPortablePdbStream(pdb);
    MetadataReader reader = provider.GetMetadataReader();
    using FileStream image = File.OpenRead(args[1]);
    using var executable = new PEReader(image);
    DebugDirectoryEntry[] codeViews = executable.ReadDebugDirectory()
        .Where(entry => entry.Type == DebugDirectoryEntryType.CodeView && entry.IsPortableCodeView).ToArray();
    if (codeViews.Length != 1 || reader.DebugMetadataHeader == null)
        throw new BadImageFormatException("Exactly one portable CodeView record and PDB identity are required.");
    CodeViewDebugDirectoryData codeView = executable.ReadCodeViewDebugDirectoryData(codeViews[0]);
    var pdbIdentity = new BlobContentId(reader.DebugMetadataHeader.Id);
    if (codeView.Guid != pdbIdentity.Guid || codeViews[0].Stamp != pdbIdentity.Stamp || codeView.Age != 1)
        throw new BadImageFormatException("The portable PDB identity does not match the assembly's CodeView record.");
    // Portable PDB custom-debug-information kinds defined by the .NET format specification.
    var sourceLinkKind = new Guid("cc110556-a091-4d38-9fec-25ab9a351a6a");
    var embeddedSourceKind = new Guid("0e8a571b-6926-466e-b4ad-8ab04611f5fe");
    JsonElement? sourceLink = null;
    foreach (CustomDebugInformationHandle handle in reader.CustomDebugInformation)
    {
        CustomDebugInformation information = reader.GetCustomDebugInformation(handle);
        if (information.Parent.Kind != HandleKind.ModuleDefinition || reader.GetGuid(information.Kind) != sourceLinkKind)
            continue;
        if (sourceLink.HasValue)
            throw new BadImageFormatException("The portable PDB contains multiple Source Link maps.");
        using JsonDocument map = JsonDocument.Parse(reader.GetBlobBytes(information.Value));
        sourceLink = map.RootElement.Clone();
    }

    var documents = new List<object>();
    foreach (DocumentHandle handle in reader.Documents.OrderBy(handle => reader.GetString(reader.GetDocument(handle).Name), StringComparer.Ordinal))
    {
        Document document = reader.GetDocument(handle);
        bool embeddedSource = reader.GetCustomDebugInformation(handle)
            .Any(information => reader.GetGuid(reader.GetCustomDebugInformation(information).Kind) == embeddedSourceKind);
        documents.Add(new
        {
            name = reader.GetString(document.Name),
            hashAlgorithm = reader.GetGuid(document.HashAlgorithm).ToString("D", CultureInfo.InvariantCulture),
            hash = Convert.ToHexString(reader.GetBlobBytes(document.Hash)),
            embeddedSource
        });
    }
    Console.Write(JsonSerializer.Serialize(new { schemaVersion = 1, sourceLink, documents },
        new JsonSerializerOptions { WriteIndented = true, NewLine = "\n" }) + "\n");
    return 0;
}

Assembly assembly = Assembly.LoadFrom(Path.GetFullPath(args[1]));
if (string.Equals(args[0], "coverage-types", StringComparison.Ordinal))
{
    using FileStream pdb = File.OpenRead(Path.ChangeExtension(args[1], ".pdb"));
    using MetadataReaderProvider provider = MetadataReaderProvider.FromPortablePdbStream(pdb);
    MetadataReader reader = provider.GetMetadataReader();
    SortedSet<string> names = new(StringComparer.Ordinal);
    foreach (MethodDebugInformationHandle handle in reader.MethodDebugInformation)
    {
        if (!reader.GetMethodDebugInformation(handle).GetSequencePoints().Any(point => !point.IsHidden))
            continue;
        MethodBase? method = assembly.ManifestModule.ResolveMethod(0x06000000 | MetadataTokens.GetRowNumber(handle));
        Type? type = method?.DeclaringType;
        while (type?.DeclaringType != null && type.IsDefined(typeof(System.Runtime.CompilerServices.CompilerGeneratedAttribute)))
            type = type.DeclaringType;
        if (type?.FullName != null)
            names.Add(type.FullName.Replace('+', '/'));
    }
    Console.Write(JsonSerializer.Serialize(names, new JsonSerializerOptions { WriteIndented = true, NewLine = "\n" }) + "\n");
    return 0;
}

List<string> lines = new();
NullabilityInfoContext nullability = new();
const BindingFlags flags = BindingFlags.Public | BindingFlags.Instance | BindingFlags.Static | BindingFlags.DeclaredOnly;
foreach (Type type in assembly.GetExportedTypes().OrderBy(type => type.FullName, StringComparer.Ordinal))
{
    string name = TypeName(type);
    lines.Add($"type {name} base={TypeName(type.BaseType)} abstract={type.IsAbstract} sealed={type.IsSealed} interfaces={string.Join(",", type.GetInterfaces().Select(TypeName).Order(StringComparer.Ordinal))}");
    foreach (ConstructorInfo constructor in type.GetConstructors(flags))
        lines.Add($"constructor {name}({Parameters(constructor.GetParameters())})");
    foreach (MethodInfo method in type.GetMethods(flags).Where(method => !method.IsSpecialName))
        lines.Add(MethodSignature(name, method));
    foreach (MethodInfo method in type.GetMethods(flags).Where(method => method.Name.StartsWith("op_", StringComparison.Ordinal)))
        lines.Add(MethodSignature(name, method));
    foreach (PropertyInfo property in type.GetProperties(flags))
        lines.Add($"property {name}.{property.Name}({Parameters(property.GetIndexParameters())}):{TypeName(property.PropertyType)} nullability={nullability.Create(property).ReadState}/{nullability.Create(property).WriteState} get={Accessor(property.GetMethod)} set={Accessor(property.SetMethod)}");
    foreach (FieldInfo field in type.GetFields(flags))
        lines.Add($"field {name}.{field.Name}:{TypeName(field.FieldType)} nullability={nullability.Create(field).ReadState}/{nullability.Create(field).WriteState} static={field.IsStatic} readonly={field.IsInitOnly} value={(field.IsLiteral ? Convert.ToString(field.GetRawConstantValue(), CultureInfo.InvariantCulture) : "-")}");
    foreach (EventInfo item in type.GetEvents(flags))
        lines.Add($"event {name}.{item.Name}:{TypeName(item.EventHandlerType)}");
}
Console.Write(string.Join("\n", lines.Order(StringComparer.Ordinal)) + "\n");
return 0;

static string TypeName(Type? type) => type?.FullName ?? type?.ToString() ?? "-";
static string Accessor(MethodInfo? method) => method?.IsPublic == true ? $"public,static={method.IsStatic},virtual={method.IsVirtual},final={method.IsFinal}" : "-";
static string Parameters(ParameterInfo[] parameters) => string.Join(",", parameters.Select(parameter =>
    $"{TypeName(parameter.ParameterType)} {parameter.Name} nullability={new NullabilityInfoContext().Create(parameter).ReadState}/{new NullabilityInfoContext().Create(parameter).WriteState} in={parameter.IsIn} out={parameter.IsOut} optional={parameter.IsOptional}" +
    (parameter.HasDefaultValue ? $" default={Convert.ToString(parameter.DefaultValue, CultureInfo.InvariantCulture) ?? "null"}" : "")));
static string MethodSignature(string name, MethodInfo method) =>
    $"method {name}.{method.Name}({Parameters(method.GetParameters())}):{TypeName(method.ReturnType)} nullability={new NullabilityInfoContext().Create(method.ReturnParameter).ReadState} static={method.IsStatic} virtual={method.IsVirtual} abstract={method.IsAbstract} final={method.IsFinal}";
