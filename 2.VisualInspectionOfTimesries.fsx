System.Environment.CurrentDirectory <- __SOURCE_DIRECTORY__

#I @".NET\PulseEngineApp\bin\Debug"
#r "Angara.table.dll"
#r "Angara.Statistics.dll"
#r "Angara.Html.dll"

open Angara
open Angara.Data
open Angara.Statistics

let ch_data = Table.Load("1.ChAreaApproximated.csv")

Html.Save "Observations.html" ch_data