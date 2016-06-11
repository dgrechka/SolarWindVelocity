let source = __SOURCE_DIRECTORY__

#I "../packages/FSharp.Formatting.2.0.2/lib/net40"
#r "FSharp.CodeFormat.dll"
#r "FSharp.Literate.dll"
open FSharp.Literate
open System.IO

let doc = Path.Combine(source, "../../MCMCFiting.md")
Literate.ProcessMarkdown(doc, template)