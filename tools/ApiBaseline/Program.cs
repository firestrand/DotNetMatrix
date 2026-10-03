using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata;
using System.Reflection.Metadata.Ecma335;
using System.Text.Json;

if (args.Length != 2 || (args[0] != "api" && args[0] != "coverage-types"))
{
    Console.Error.WriteLine("Usage: ApiBaseline api|coverage-types ASSEMBLY_PATH");
    return 2;
}

Assembly assembly = Assembly.LoadFrom(Path.GetFullPath(args[1]));
if (args[0] == "coverage-types")
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
    Console.WriteLine(JsonSerializer.Serialize(names, new JsonSerializerOptions { WriteIndented = true }));
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
Console.WriteLine(string.Join(Environment.NewLine, lines.Order(StringComparer.Ordinal)));
return 0;

static string TypeName(Type? type) => type?.FullName ?? type?.ToString() ?? "-";
static string Accessor(MethodInfo? method) => method?.IsPublic == true ? $"public,static={method.IsStatic},virtual={method.IsVirtual},final={method.IsFinal}" : "-";
static string Parameters(ParameterInfo[] parameters) => string.Join(",", parameters.Select(parameter =>
    $"{TypeName(parameter.ParameterType)} {parameter.Name} nullability={new NullabilityInfoContext().Create(parameter).ReadState}/{new NullabilityInfoContext().Create(parameter).WriteState} in={parameter.IsIn} out={parameter.IsOut} optional={parameter.IsOptional}" +
    (parameter.HasDefaultValue ? $" default={Convert.ToString(parameter.DefaultValue, CultureInfo.InvariantCulture) ?? "null"}" : "")));
static string MethodSignature(string name, MethodInfo method) =>
    $"method {name}.{method.Name}({Parameters(method.GetParameters())}):{TypeName(method.ReturnType)} nullability={new NullabilityInfoContext().Create(method.ReturnParameter).ReadState} static={method.IsStatic} virtual={method.IsVirtual} abstract={method.IsAbstract} final={method.IsFinal}";
