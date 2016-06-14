﻿module app

open BurstEngine

open Angara
open Angara.Data
open Angara.Statistics
open Angara.Filzbach

open MathNet.Numerics.LinearAlgebra

//tick is 1 minute real time
let timeToTick t =
    int(round(t*60.0))

let a = 1000.0 * 60.0 / 10.0**8.0
let velocityToModel v = v * a
let toKmPerS v = v / a

let max_hole= 20.0
let max_velocity = 1.0

let getWind (predictorsT: Table) v_max v_p d_max d_p = 
    //velocity slope definition is transforemed
    //parameters
    //max_v (0.0 - 1.0) defained as "wind velocity generated by hole `max_hole`"
    //base_v_p (0.0 - 1.0)defined as "max_v * base_v_p = wind velocity generate by hole near zero area"

    let generatePulse (t:float) area =
        let emergence_tick = timeToTick(t)
        let area = if area<0.0 then 0.0 else area
        if area = 0.0 then
            None
        else
            //max area is considered to be ~ 20
            //should give velocity of 1.0
            //thus
            let k = v_max * (1.0 - v_p) / max_hole
            let b = v_max * v_p
            let velocity = k * area + b                    
            #if DEBUG
            assert(abs(k* max_hole + b - v_max)<1e-14)
            assert(abs(k*0.0 + b - v_max*v_p)<1e-14)
            #endif
            let k = d_max * (1.0 - d_p) / max_hole
            let b = d_max * d_p
            let density = k * area + b                    
            #if DEBUG
            assert(abs(k* max_hole + b - d_max)<1e-14)
            assert(abs(k*0.0 + b - d_max*d_p)<1e-14)
            #endif            
            let burst = {
                EmergenceTime=emergence_tick;
                Velocity=velocity;
                Density= density
            }
            Some(burst)
    Table.Map ["ts";"Relative_CH_CorrectSphereArea_j"] generatePulse predictorsT |> Seq.choose (fun b -> b)  |> List.ofSeq
  

let mutable bestLglk = System.Double.NegativeInfinity

let logLikelihood (predictorsT: Table) (observationsT: Table) (p: Parameters) =                
    //Distance to Earth
    //varies between 0.9832898912 and 1.0167103335 AU. 
    //which is 1.47098074 to 1.52097701 ×10^8 km
    
    //1 model space step is 1.0x10^8m which is approx 10 Earth diameter
    //thus 1 AU is ~ 1500 model space steps

    let D = p.GetValue "D" //distance to Earth
    let Sigma = p.GetValue "Sigma" //observation noise    
    let v_max = p.GetValue "V_max" //wind velocity generated by constant max_hole
    let v_p = p.GetValue "V_p"
    let d_max = p.GetValue "D_max" //wind density generated by constant max_hole
    let d_p = p.GetValue "D_p"
    let d_b = p.GetValue "D_bg" //density of background wind
    let v_b = p.GetValue "V_bg" //velocity of background wind

    let wind = getWind predictorsT v_max v_p d_max d_p |> expandWind    
    
    let data = Table.Map ["ts";"velocity_mean"] (fun (t:float) (v:float) -> (t,v)) observationsT    

    let current_lglk p =
        let t,v = p
        let v =  velocityToModel v //velocity to model units (1.0 = 1.0x10^8m/min)
        let t = timeToTick t //hours to model units (minutes)
        let burstsAtEarth = locationState wind t (int(round(D)))
        let windAtEart = windAvg burstsAtEarth
        let v_avg_Earth,d_avg_Earth = windAtEart
        let mu = (v_b*d_b+v_avg_Earth*d_avg_Earth)/(d_b+d_avg_Earth)
        if System.Double.IsNaN mu then
            log improbable
        else
            let distribution = Statistics.Normal(mu,Sigma)
            let lglk = log_pdf distribution v
            //printfn "mu:%g\tv:%g\tsigma:%g\tlglk:%g" mu v Sigma lglk
            lglk        
    let res =
        data |> Seq.toArray
//             |> Array.Parallel.map current_lglk
             |> Array.map current_lglk
             |> Array.sum         
    if res>bestLglk then
        printfn "\nlglk improvement:%g\tVm:%g\tVp:%g\tDm:%g\tDp:%g\tVbg:%g\tDbg:%g\tSig:%g\tD:%g"
            res v_max v_p d_max d_p v_b d_b Sigma D
        bestLglk <- res
    else
        printf "."
    res
        

[<EntryPoint>]
let main argv =     
    //config
    let doEstimate = false
    let thinnObs = false
    let seed : uint32 ref = ref 0u
    if argv.Length>0
        then System.UInt32.TryParse(argv.[0],seed) |> ignore
    let rng = new MT19937(!seed)


    //action
    printfn "Seed is %d" !seed

    let ch_data = Table.Load @"..\..\..\..\ResultData\CH_features_cleaned_2015.csv"
    let obs = Table.Load @"..\..\..\..\ResultData\ACE_EPAM_SW_2015.csv"    
    
    let mutable counter = 0
    let obs =
        if thinnObs then
            Table.Filter ["ts"]  (fun (dummy:float) ->            
            counter <- counter+1
            counter%5=0            
            ) obs
        else obs

    let predictorsT = Table.Filter ["ts"] (fun t -> t > 40.0 && t < 2140.0) ch_data
    let observationsT = Table.Filter ["ts"] (fun t -> t > 150.0 && t < 2140.0) obs

    printfn "Predictor values count %d" predictorsT.RowsCount
    printfn "Observation values count %d" observationsT.RowsCount    
    
    let eps = System.Double.Epsilon

    let estimate predictorsT observationsT =
        let parameters = Parameters.Empty
                            .Add("V_max",[|0.8|],eps,1.0,isLog=true,delay=0)
                            .Add("V_p",[|0.5|],eps,1.0,isLog=true,delay=0)
                            .Add("D_max",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("D_p",[|0.5|],eps,1.0,isLog=true,delay=0)
                            .Add("V_bg",[|0.1|],eps,0.5,isLog=true,delay=0)
                            .Add("D_bg",[|0.5|],eps,1.0,isLog=false,delay=0)
                            //.Add("D_bg",[|0.5|])
                            .Add("D",[|1500.0|],1200.0,1800.0,isLog=false,delay=0)
                            .Add("Sigma",[|0.01|],eps, 0.5,isLog=true,delay=0)                            
        Sampler.runmcmc(parameters, logLikelihood predictorsT observationsT, 10000, 5000,rng=rng)
    if doEstimate then
        let res = estimate predictorsT observationsT    
        Sampler.print res
//        let posterior  =
//            res.samples 
//            |> Seq.mapi (fun i {values=sample} -> Column.Create(res.sampler.Parameters.GetName(i),sample))
//            |> Table.OfColumns
//        printfn "done"
        //Table.Save(posterior, sprintf "posterior_%d.csv" !seed)        
    else        
        //simulation    
        let v_m = 0.550798
        let v_p = 0.219125
        let d_m = 0.883239
        let d_p = 7.6e-49
        let v_b = 0.269306
        let d_b = 0.469225
        let D = 1587.39
        let Sigma = 0.0430126

        let testWind = getWind predictorsT v_m v_p d_m d_p

        let testWind = expandWind testWind

        let time_start = timeToTick 150.0        
        let time_stop = timeToTick 200.0        
        let space_start = 1400
        let space_end = 1450

        let simulation = simulate testWind time_start time_stop space_start space_end

        //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
        let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
        ds.IsAutocommitEnabled <- false
        //let windVar = ds.AddVariable<float>("windDens",[|"x";"t"|])
        let windSpVar = ds.AddVariable<float>("avgWindSpeed",[|"x";"t"|])
        //let windSpMaxVar = ds.AddVariable<float>("maxWindSpeed",[|"x";"t"|])
        let timeAxis = ds.AddVariable<int>("t",[|"t"|])
        let spaceAxis = ds.AddVariable<int>("x",[|"x"|])        

        let timeObsAxis = ds.AddVariable<int>("obs_time",[|"t_obs"|])
        let obsAxis = ds.AddVariable<float>("obs",[|"t_obs"|])
        let predAxis = ds.AddVariable<float>("pred",[|"t_obs"|])
        let predMeanAxis = ds.AddVariable<float>("predMean",[|"t_obs"|])

        let dummy =
            Table.MapToColumn "pred" ["ts";"velocity_mean"]  (fun (t:float) (v:float) ->
                let burstsAtEarth = locationState testWind (timeToTick t) (int(round(D)))
                let windAtEart = windAvg burstsAtEarth
                let v_avg_Earth,d_avg_Earth = windAtEart
                let mu = (v_b*d_b+v_avg_Earth*d_avg_Earth)/(d_b+d_avg_Earth)
                        
                let n = Statistics.Normal(mu,Sigma)
                let pred = Statistics.draw rng n
                timeObsAxis.Append([|timeToTick t|]);
                obsAxis.Append([|v|]);
                predAxis.Append([|toKmPerS pred|]);
                predMeanAxis.Append([|toKmPerS mu|]);
                //printfn "%f %f %f %f" t v pred mu
                pred
                ) observationsT

        printfn "simulated for %d time moments" (dummy.RowsCount)

        List.iteri (fun i sample_speed ->
            let t = time_start+i
            let sample_speed_km_s = Array.map toKmPerS sample_speed
            windSpVar.Append(sample_speed_km_s,"t")            
            timeAxis.Append([|t|]);
            )
            simulation
        let spaceAxisData = Array.init (space_end-space_start+1) (fun i -> space_start+i)
        spaceAxis.Append(spaceAxisData)

        ds.Commit()        
    
    0 // return an integer exit code
