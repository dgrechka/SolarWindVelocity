module app

open BurstEngine
open Angara.Data

let loadPulses (csvFile:string) =
    Table.Load(csvFile) |>
    Table.Map ["t";"density";"velocity"] (fun (t:float) (d:float) (v:float) ->
        {
            EmergenceTime=int(t);
            Density=d;
            Velocity=v            
        }
    ) |> List.ofSeq

[<EntryPoint>]
let main argv = 
    Angara.Base.Init()
    
    let testWind = loadPulses @"..\..\..\..\TestData\6pulses.csv"
    
    let testWind = expandWind testWind

    let space_start = 10
    let space_end = 120
    let time_start = 0
    let time_end = 120
    
    let simulation = simulate testWind time_start time_end 3 space_start space_end

    //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
    let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
    ds.IsAutocommitEnabled <- false    
    let windSpVar = ds.AddVariable<float>("avgWindSpeed",[|"x";"t"|])    
    let timeAxis = ds.AddVariable<int>("t",[|"t"|])
    let spaceAxis = ds.AddVariable<int>("x",[|"x"|])    

    List.iteri (fun i sampls_speed ->        
        let t = i        
        windSpVar.Append(sampls_speed,"t")        
        timeAxis.Append([|t|]);
        )
        simulation
    let spaceAxisData = Array.init (space_end - space_start + 1) (fun i -> space_start+i)
    spaceAxis.Append(spaceAxisData)

    ds.Commit()        
    0 // return an integer exit code
