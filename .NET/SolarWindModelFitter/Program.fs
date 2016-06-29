module app

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
      
let get2ndPolyCoefs x0 x1 x2 y0 y1 y2 = 
    let m = matrix [[ x0*x0; x0; 1.0]
                    [ x1*x1; x1; 1.0]
                    [ x2*x2; x2; 1.0]]
    let r_part = vector [y0; y1; y2]
    let solution = m.Solve r_part
    assert(Vector.length solution = 3)
    (solution.[0],solution.[1],solution.[2]) 

let getWind (predictorsT: Table) v1p v2p v3p d1p d2p d3p vap dap= 
    //velocity curve definition is transforemed
    //all parameters are portions:

    //velocity(area) function crosses these 3 points:
    //v3 = velocity(max_hole) = max_velocity * v3p
    //v2 = velocity(ap) = v3 * v2p
    //v1 = velocity(0.0) = v2 * v1p

    //denity(area) function crosses these 3 points:
    //d3 = density(max_hole) =  d3p
    //d2 = density(ap) = d3 * d2p
    //d1 = density(0.0) = d2 * d1p

    let v2 = v3p*v2p
    let d2 = d3p*d2p
    let va,vb,vc = get2ndPolyCoefs 0.0 (max_hole*vap) max_hole (v2*v1p) v2 v3p
    let da,db,dc = get2ndPolyCoefs 0.0 (max_hole*dap) max_hole (d2*d1p) d2 d3p
            

    let generatePulse (t:float) area =
        let emergence_tick = timeToTick(t)
        let area = if area<0.0 then 0.0 else area
        if area = 0.0 then
            None
        else
            let areaSq = area*area
            let burst = {
                EmergenceTime=emergence_tick;
                Velocity= va*areaSq + vb*area + vc;
                Density= da*areaSq + db*area + dc
            }
            Some(burst)
    Table.Map ["ts";"Relative_CH_CorrectSphereArea_j"] generatePulse predictorsT |> Seq.choose (fun b -> b)  |> List.ofSeq
  

let mutable bestLglk = System.Double.NegativeInfinity

let mutable iteration_counter = 0

let logLikelihood (predictorsT: Table) (observationsT: Table) (p: Parameters) =                
    //Distance to Earth
    //varies between 0.9832898912 and 1.0167103335 AU. 
    //which is 1.47098074 to 1.52097701 ×10^8 km
    
    //1 model space step is 1.0x10^8m which is approx 10 Earth diameter
    //thus 1 AU is ~ 1500 model space steps

    let D = p.GetValue "D" //distance to Earth
    let Sigma = p.GetValue "Sigma" //observation noise    
    let v3p = p.GetValue "V_p3"
    let v2p = p.GetValue "V_p2"
    let v1p = p.GetValue "V_p1"
    let d3p = p.GetValue "D_p3"
    let d2p = p.GetValue "D_p2"
    let d1p = p.GetValue "D_p1"
    let vap = p.GetValue "Area_vp"
    let dap = p.GetValue "Area_dp"
    let d_b = p.GetValue "D_bg" //density of background wind
    let v_b = p.GetValue "V_bg" //velocity of background wind

    let wind = getWind predictorsT v1p v2p v3p d1p d2p d3p vap dap |> expandWind    
    
    let data = Table.Map ["ts";"velocity_mean"] (fun (t:float) (v:float) -> (t,v)) observationsT    

    let current_lglk p wind =
        let t,v = p
        let v =  velocityToModel v //velocity to model units (1.0 = 1.0x10^8m/min)
        let t = timeToTick t //hours to model units (minutes)
        let burstsAtEarth,optimizedWind = locationStateOptimized wind t (int(round(D)))
        let windAtEart = windAvg burstsAtEarth
        let v_avg_Earth,d_avg_Earth = windAtEart
        let mu = (v_b*d_b+v_avg_Earth*d_avg_Earth)/(d_b+d_avg_Earth)
        let result =
            if System.Double.IsNaN mu then
                log improbable
            else
                let distribution = Statistics.Normal(mu,Sigma)
                //let varSq = Sigma*Sigma*Sigma*Sigma
                //let alpha_shape = mu*mu/varSq
                //let beta_rate = mu/varSq
                //let distribution = Statistics.Gamma(alpha_shape,beta_rate)                         
                let lglk = log_pdf distribution v
                //printfn "mu:%g\tv:%g\tsigma:%g\tlglk:%g" mu v Sigma lglk
                lglk        
        result,optimizedWind
    
    let folder state observation =
        let lglk_sum,wind = state
        let curLgLk,optimizedWind = current_lglk observation wind
        //printfn "was %d now %d" (List.length wind) (List.length optimizedWind)
        (lglk_sum+curLgLk),optimizedWind

    let folded =
        data |> Seq.toArray |>
        Array.fold folder (0.0,wind)
    
    let res = fst folded
                      
    iteration_counter <- iteration_counter + 1
    if res>bestLglk then        
        printfn "\n----------"
        printfn "\niteration:%d\tlog-likelihood:%g\timprovement:%g" iteration_counter res (res-bestLglk)
        printfn "Vfunc nodes:\t(%g;%g) (%g;%g) (%g;%g)" 0.0 (v3p*v2p*v1p) (vap*max_hole) (v3p*v2p) max_hole v3p
        printfn "Dfunc nodes:\t(%g;%g) (%g;%g) (%g;%g)"0.0 (d3p*d2p*d1p) (dap*max_hole) (d3p*d2p) max_hole d3p
        printfn "background wind\tV:%g\tD:%g" v_b d_b
        printfn "Noise sigma:\t%g" Sigma
        printfn "Earth distance:\t%f" D
        System.Console.Beep(3500,500)
        bestLglk <- res
    else
        printf "."
        System.Console.Out.Flush()
        System.Console.Beep(2000,100)
        //printfn "lglk: %g" res

    res    

[<EntryPoint>]
let main argv =     
    //tests
    let (a,b,c) = get2ndPolyCoefs -1.0 0.0 1.0 2.0 1.0 2.0
    assert(a=1.0)
    assert(b=0.0)
    assert(c=1.0)

    let (a,b,c) = get2ndPolyCoefs 0.0 1.0 2.0 2.0 1.0 6.0
    assert(a=3.0)
    assert(b= -4.0)
    assert(c=2.0)
    printfn "Basic tests passed"

    //config
    let doEstimate = true
    let thinnObs = false
    let seed : uint32 ref = ref 0u
    if argv.Length>0
        then System.UInt32.TryParse(argv.[0],seed) |> ignore
    let rng = new MT19937(!seed)


    //action
    printfn "Seed is %d" !seed

    let ch_data = Table.Load @"CH_features_cleaned_2015.csv"
    let obs = Table.Load @"ACE_EPAM_SW_2015.csv"    
    
    let mutable counter = 0
    let obs =
        if thinnObs then
            Table.Filter ["ts"]  (fun (dummy:float) ->            
            counter <- counter+1
            counter%5=0            
            ) obs
        else obs

    let predictorsT = Table.Filter ["ts"] (fun t -> t > 40.0 && t < 8000.0) ch_data
    let observationsT = Table.Filter ["ts"] (fun t -> t > 150.0 && t < 8000.0) obs

    printfn "Predictor values count %d" predictorsT.RowsCount
    printfn "Observation values count %d" observationsT.RowsCount    
    
    let eps = System.Double.Epsilon

    let estimate predictorsT observationsT =
        let parameters = Parameters.Empty
                            .Add("V_p3",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("V_p2",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("Area_vp",[|0.5|],0.1,1.0,isLog=false,delay=0)
                            .Add("V_p1",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("D_p3",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("D_p2",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("Area_dp",[|0.5|],0.1,1.0,isLog=false,delay=0)
                            .Add("D_p1",[|0.8|],eps,1.0,isLog=false,delay=0)                            
                            .Add("V_bg",[|0.1|],eps,0.5,isLog=true,delay=0)
                            .Add("D_bg",[|0.5|],eps,10.0,isLog=true,delay=0)
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
        0     
    else        
        //simulation    
//        let v_m = 0.470203
//        let v_p = 0.558267
//        let d_m = 0.519722
//        let d_p = 1.25e-245
//        let v_b = 0.224688
//        let d_b = 0.675671
//        let D = 1795.86
//        let Sigma = 0.0527357
//
//        let testWind = getWind predictorsT v_m v_p d_m d_p
//
//        let testWind = expandWind testWind
//
//        let time_start = timeToTick 1500.0        
//        let time_stop = timeToTick 2000.0        
//        let time_by = 100
//        let space_start = 0
//        let space_end = 1587
//
//        let simulation = simulate testWind time_start time_stop time_by space_start space_end
//
//        //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
//        let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
//        ds.IsAutocommitEnabled <- false
//        //let windVar = ds.AddVariable<float>("windDens",[|"x";"t"|])
//        let windSpVar = ds.AddVariable<float>("avgWindSpeed",[|"x";"t"|])
//        //let windSpMaxVar = ds.AddVariable<float>("maxWindSpeed",[|"x";"t"|])
//        let timeAxis = ds.AddVariable<int>("t",[|"t"|])
//        let spaceAxis = ds.AddVariable<int>("x",[|"x"|])        
//
//        let timeObsAxis = ds.AddVariable<int>("obs_time",[|"t_obs"|])
//        let obsAxis = ds.AddVariable<float>("obs",[|"t_obs"|])
//        let predAxis = ds.AddVariable<float>("pred",[|"t_obs"|])
//        let predMeanAxis = ds.AddVariable<float>("predMean",[|"t_obs"|])
//
//        let dummy =
//            Table.MapToColumn "pred" ["ts";"velocity_mean"]  (fun (t:float) (v:float) ->
//                let burstsAtEarth = locationState testWind (timeToTick t) (int(round(D)))
//                let windAtEart = windAvg burstsAtEarth
//                let v_avg_Earth,d_avg_Earth = windAtEart
//                let mu = (v_b*d_b+v_avg_Earth*d_avg_Earth)/(d_b+d_avg_Earth)
//                        
//                let n = Statistics.Normal(mu,Sigma)
//                let pred = Statistics.draw rng n
//                timeObsAxis.Append([|timeToTick t|]);
//                obsAxis.Append([|v|]);
//                predAxis.Append([|toKmPerS pred|]);
//                predMeanAxis.Append([|toKmPerS mu|]);
//                //printfn "%f %f %f %f" t v pred mu
//                pred
//                ) observationsT
//
//        printfn "simulated for %d time moments" (dummy.RowsCount)
//
//        List.iteri (fun i sample_speed ->
//            let t = time_start+i
//            let sample_speed_km_s = Array.map toKmPerS sample_speed
//            windSpVar.Append(sample_speed_km_s,"t")            
//            timeAxis.Append([|t|]);
//            )
//            simulation
//        let spaceAxisData = Array.init (space_end-space_start+1) (fun i -> space_start+i)
//        spaceAxis.Append(spaceAxisData)
//
//        ds.Commit()        
    
    0 // return an integer exit code
