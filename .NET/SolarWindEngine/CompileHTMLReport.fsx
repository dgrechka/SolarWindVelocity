let source = __SOURCE_DIRECTORY__

#I "../packages/FSharp.Formatting.2.14.4/lib/net40"
#I "bin\\Debug"
#r "FSharp.CodeFormat.dll"
#r "FSharp.Literate.dll"
open FSharp.Literate
open System.IO

let projInfo =
  [ "page-description", "F# Literate programming"
    "project-author", "Tomas Petricek"
    "github-link", "https://github.com/tpetricek/FSharp.Formatting"
    "project-name", "F# Formatting" ]

let template = Path.Combine(source, "../packages/FSharp.Formatting.2.14.4/templates/template.cshtml")
let script = Path.Combine(source, "PulseEngine.fsx")
Literate.ProcessScriptFile(script,template, replacements = projInfo)