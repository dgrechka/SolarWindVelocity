module app

open PulseEngine
open Angara
open Angara.Data
open Angara.Statistics
open Angara.Filzbach

let getWind (predictorsT: Table) S_0 V_k P_w V_b = 
    let generatePulse t area = 
        if area < S_0 then
            None
        else
            let intercept = - V_k*S_0
            let A = intercept + V_k*area
            let velocity = V_b+A
            assert(A >= 0.0)
            let pulse = {
                Power=1.0;
                EmergenceTime=t;
                Velocity=velocity;
                A=A;
                Kernel=Kernels.GaussianExt P_w
            }
            Some(pulse)
    Table.Map ["t";"CHarea"] generatePulse predictorsT |> Seq.choose (fun a -> a) |> List.ofSeq

let logLikelihood (predictorsT: Table) (observationsT: Table) (p: Parameters) =        
    let V_k = p.GetValue "V_k" //velocity growth slope
    let S_0 = p.GetValue "S_0" //minimal CH area that enables slope
    let V_b = p.GetValue "V_b" //"background" velocity

    //Distance to Earth
    //varies between varies between 0.9832898912 and 1.0167103335 AU. 
    //the same as 1.47098074 to 1.52097701 ×10^8 km
    let D = p.GetValue "D" //distance to Earth
    let Sigma = p.GetValue "Sigma" //observation noise
    let P_w = p.GetValue "P_w" //Pulse width. around pulse distance traveled in half an hour?


    let wind = getWind predictorsT S_0 V_k P_w V_b
        
    let lglks: float seq = 
        let current_lglk t v =
            let avgVelocity = WindAvgAForTime wind t
            let predAtEarth = avgVelocity D
            let mu = V_b + predAtEarth //base level plus pulses
            let distribution = Statistics.Normal(mu,Sigma)
            let lglk = log_pdf distribution v
            lglk
        Table.Map ["t";"velocity"] current_lglk observationsT
    lglks |> Seq.sum       

[<EntryPoint>]
let main argv =         
    let ch_data = Table.Load @"..\..\..\..\1.ChAreaApproximated.csv"
    let obs = Table.Load @"..\..\..\..\SampleData\ace_swepam_2015_1h.csv"
    let obs2 = Table.MapToColumn ["hours since 01.01.2015  1:00:00"] "t" (fun t -> t+1.0) obs

    let predictorsT = Table.Filter ["t"] (fun t -> t > 40.0 ) ch_data
    let observationsT = Table.Filter ["t"] (fun t -> t > 40.0+168.0 ) obs2 //first 168 hours is burnin to enable first pulse reach the Earth

    let time_start = 0.0
    let time_step = 0.2
    let time_stop = 700.0
    let space_step = 0.5e+6
    let space_start = 0.0 
    let space_end = 1.5e+8        

    
    //simulation
    
    let testFastWind = getWind predictorsT 0.0 (40.0*60.0*60.0) (300.0*60.0*120.0) (300.0*60.0*60.0)    

    let simulation = simulate testFastWind time_start time_stop time_step space_start space_end space_step

    //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
    let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
    ds.IsAutocommitEnabled <- false
    let windVar = ds.AddVariable<float>("windDens",[|"x";"t"|])
    let windSpVar = ds.AddVariable<float>("windSpeed",[|"x";"t"|])
    let timeAxis = ds.AddVariable<float>("t",[|"t"|])
    let spaceAxis = ds.AddVariable<float>("x",[|"x"|])        

    List.iteri (fun i sim_step_data ->
        let sample_vals,sample_speed = sim_step_data
        let sample_speed_km_s = List.map (fun a -> a/60.0/60.0 + 300.0) sample_speed
        let t = time_start+float(i)*time_step
        windVar.Append(List.toArray sample_vals,"t")
        windSpVar.Append(List.toArray sample_speed_km_s,"t")
        timeAxis.Append([|t|]);
        )
        simulation
    let spaceAxisData = Array.init (int((space_end-space_start)/space_step)) (fun i -> space_start+float(i)*space_step)
    spaceAxis.Append(spaceAxisData)

    ds.Commit()        
    
    0 // return an integer exit code
