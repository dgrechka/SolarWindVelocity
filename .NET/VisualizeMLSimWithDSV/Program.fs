// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.

open Microsoft.Research.Science.Data

[<EntryPoint>]
let main argv = 
    use ds = DataSet.Open @"msds:csv?file=MLSimulation.csv&openMode=readOnly&inferDims=true"
    use outDS = DataSet.Open @"msds:nc?file=out.nc&openMode=create"
    outDS.IsAutocommitEnabled <- false;
    let futureVar = outDS.AddVariable<float>("upcoming",[|"t";"fut_space"|]);
    let pastVar = outDS.AddVariable<float>("past",[|"t";"past_space"|]);
    let currentVar = outDS.AddVariable<float>("current",[|"t";"current_space"|]);
    let futurePredVar = outDS.AddVariable<float>("pred_m",[|"t";"pred_space"|]);
    let futureUpperVar = outDS.AddVariable<float>("pred_u",[|"t";"pred_space"|]);
    let futureLowerVar = outDS.AddVariable<float>("pred_l",[|"t";"pred_space"|]);
    let timeVar = outDS.AddVariable<float>("time",[|"t"|]);
    let N = ds.Dimensions.[0].Length

    let fut_space_var = outDS.AddVariable<float>("fut_space_var",[|"fut_space"|]);
    fut_space_var.Append(Array.init 47 (fun i -> float (i+1)));
    let past_space_var = outDS.AddVariable<float>("past_space_var",[|"past_space"|]);
    past_space_var.Append(Array.init 143 (fun i -> float(-144+i)));
    let current_space_var = outDS.AddVariable<float>("current_space_var",[|"current_space"|]);
    current_space_var.Append([|0.0|]);
    let pred_space_var = outDS.AddVariable<float>("pred_space_var",[|"pred_space"|]);
    pred_space_var.Append(Array.init 16 (fun i -> ((float i)+1.0)*3.0));
    fut_space_var.Metadata.["Units"] <- "hours";
    current_space_var.Metadata.["Units"] <- "hours";
    past_space_var.Metadata.["Units"] <- "hours";
    pred_space_var.Metadata.["Units"] <- "hours";

    let vel = ds.Variables.["vel"].GetData();
    let time = ds.Variables.["ts"].GetData();
    let future = Array.init 47 (fun i -> let varName = sprintf "vel.+%d" (i+1) in ds.Variables.[varName].GetData())
    let past = Array.init 143 (fun i -> let varName = sprintf "vel.-%d" (i+1) in ds.Variables.[varName].GetData())

    let pred_m = Array.init 16 (fun i -> let varName = sprintf "vel.p.+%d" (int(((float i)+1.0)*3.0)) in ds.Variables.[varName].GetData())
    let pred_sd = Array.init 16 (fun i -> let varName = sprintf "vel.p.sd.+%d" (int(((float i)+1.0)*3.0)) in ds.Variables.[varName].GetData())    

    for i in [0..N-1] do
        currentVar.Append([|vel.GetValue(i) :?> float|],0)
        timeVar.Append([|time.GetValue(i) :?> float|],0)
        let pastArray = Array.init 143 (fun j -> past.[j].GetValue(i) :?> float) |> Array.rev
        pastVar.Append(pastArray,0)
        let futureArray = Array.init 47 (fun j -> future.[j].GetValue(i) :?> float)
        futureVar.Append(futureArray,0)
        let predMArray = Array.init 16 (fun j -> pred_m.[j].GetValue(i) :?> float)
        let predUArray = Array.map2 (fun m (sd_ar:System.Array) -> m+(sd_ar.GetValue(i) :?> float)) predMArray pred_sd
        let predLArray = Array.map2 (fun m (sd_ar:System.Array) -> m-(sd_ar.GetValue(i) :?> float)) predMArray pred_sd
        futurePredVar.Append(predMArray,0)
        futureUpperVar.Append(predUArray,0)
        futureLowerVar.Append(predLArray,0)
    outDS.Commit()
    0 // return an integer exit code
