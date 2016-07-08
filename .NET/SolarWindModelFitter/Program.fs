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

//2nd order polynomial coeffecients, so the slope at (x2;y2) is 0.0
let getMaxSlope2ndPolyCoefs y0 y2 x2 = 
    assert(y0 <= y2)
    let m = matrix [[ 0.0; 0.0; 1.0]
                    [ x2*x2; x2; 1.0]
                    [ 2.0*x2; 1.0; 0.0]]
    let r_part = vector [y0; y2; 0.0]
    let solution = m.Solve r_part

    //asserting initial requrements
    assert(Vector.length solution = 3)
    let treshold = 1e-10
    assert(abs(solution.[2]-y0)<treshold)
    assert(solution.[1] > 0.0)
    assert(abs(solution.[0]*x2*x2 + solution.[1]*x2+solution.[2]-y2)<treshold)
    assert(abs(2.0*solution.[0]*x2 + solution.[1])<treshold)

    (solution.[0],solution.[1],solution.[2]) 

let get2ndPolyCoefs y0 y2 x2 k = 
    let m = matrix [[ 0.0; 0.0; 1.0]
                    [ x2*x2; x2; 1.0]
                    [ 0.0; 1.0; 0.0]]
    let r_part = vector [y0; y2; k]
    let solution = m.Solve r_part

    //asserting initial requrements
    assert(Vector.length solution = 3)
    assert(solution.[2]=y0)
    assert(solution.[0]*x2*x2 + solution.[1]*x2+solution.[2]=y2)
    assert(solution.[1]=k)

    (solution.[0],solution.[1],solution.[2]) 


let getWind (predictorsT: Table) v0p v2p v_kp d0p d2p d_kp = 
    //velocity curve definition is transforemed
    //all parameters are portions:

    //velocity(area) function crosses these 2 points:
    //v2 = velocity(max_hole) = max_velocity * v3p    
    //v0 = velocity(0.0) = v2 * v0p
    //and
    //v_k = v_kp*velocity_max_slope(0.0)

    let v2 = v2p
    let v0 = v2*v0p

    let _,vk_max,_ = getMaxSlope2ndPolyCoefs v0 v2 max_hole
    let k = vk_max*v_kp
    let va,vb,vc = get2ndPolyCoefs v0 v2 max_hole k
    
    //denity(area) function crosses these 2 points:
    //d3 = density(max_hole) =  d3p
    //d2 = density(ap) = d3 * d0p
    //and
    //d_k = d_kp*density_max_slope(0.0)

    let d2 = d2p
    let d0 = d2*d0p

    let _,dk_max,_ = getMaxSlope2ndPolyCoefs d0 d2 max_hole
    let k = dk_max*d_kp
    let da,db,dc = get2ndPolyCoefs d0 d2 max_hole k

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

let logLikelihood (predictorsT: Table) (observationsT: Table) v0p v2p v_kp d0p d2p d_kp v_b d_b D Sigma =
    //Distance to Earth
    //varies between 0.9832898912 and 1.0167103335 AU. 
    //which is 1.47098074 to 1.52097701 ×10^8 km
    
    //1 model space step is 1.0x10^8m which is approx 10 Earth diameter
    //thus 1 AU is ~ 1500 model space steps

    
    let wind = getWind predictorsT v0p v2p v_kp d0p d2p d_kp |> expandWind    
    
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
        printfn "\niteration:%d log-likelihood:%g improvement:%g" iteration_counter res (res-bestLglk)
        printfn "Vfunc nodes: (%g;%g) (%g;%g) slope: %g" 0.0 (v2p*v0p) max_hole v2p v_kp
        printfn "Dfunc nodes: (%g;%g) (%g;%g) slope: %g" 0.0 (d2p*d0p) max_hole d2p d_kp
        printfn "background wind V:%g D:%g" v_b d_b
        printfn "Noise sigma: %g" Sigma
        printfn "Earth distance: %f" D
        System.Console.Beep(3500,500)
        bestLglk <- res
    else
        printf "."
        System.Console.Out.Flush()
        System.Console.Beep(2000,100)
        //printfn "lglk: %g" res
    res

let logLikelihoodP (predictorsT: Table) (observationsT: Table) (p: Parameters) =
    let D = p.GetValue "D" //distance to Earth
    let Sigma = p.GetValue "Sigma" //observation noise        
    let v2p = p.GetValue "V_p2"
    let v0p = p.GetValue "V_p0"    
    let d2p = p.GetValue "D_p2"
    let d0p = p.GetValue "D_p0"
    let v_kp = p.GetValue "V_kp"
    let d_kp = p.GetValue "D_kp"
    let d_b = p.GetValue "D_bg" //density of background wind
    let v_b = p.GetValue "V_bg" //velocity of background wind
    logLikelihood predictorsT observationsT v0p v2p v_kp d0p d2p d_kp v_b d_b D Sigma

[<EntryPoint>]
let main argv =     
    //config
    let doEstimate = true
    let thinnObs = false
    let checkLglkForSimulation = true
    let seed : uint32 ref = ref 0u
    if argv.Length>0
        then System.UInt32.TryParse(argv.[0],seed) |> ignore
    let rng = new MT19937(!seed)


    //action
    printfn "Random seed is %d" !seed

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
                            .Add("V_p2",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("V_p0",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("V_kp",[|0.5|],eps,1.0,isLog=false,delay=0)                            
                            .Add("D_p2",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("D_p0",[|0.8|],eps,1.0,isLog=false,delay=0)
                            .Add("D_kp",[|0.5|],eps,1.0,isLog=false,delay=0)                            
                            .Add("V_bg",[|0.1|],eps,0.5,isLog=true,delay=0)
                            .Add("D_bg",[|0.5|],eps,10.0,isLog=true,delay=0)
                            .Add("D",[|1500.0|],1200.0,1800.0,isLog=false,delay=0)
                            .Add("Sigma",[|0.01|],eps, 0.5,isLog=true,delay=0)                            
        Sampler.runmcmc(parameters, logLikelihoodP predictorsT observationsT, 10000, 5000,rng=rng)
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
        ()          
        //simulation            
        (*
        let v0 = 0.228958        
        let v2 = 0.349027
        let v_kp = 0.5
        let v_bg = 3.65482e-119
        let d0 = 0.13812        
        let d2 = 0.890103        
        let d_kp = 0.5
        let d_bg = 1.15844e-123
        let D = 1225.551240
        let Sigma = 0.0497358
        let reference_lglk = 11869.9

        //converting back parameters
        let v_p2 = v2/max_velocity        
        let v_p0 = v0/v2        
        let d_p2 = d2        
        let d_p0 = d0/d2
        
        //let _,vkm,_ = getMaxSlope2ndPolyCoefs v0 v2 max_hole
                
        if checkLglkForSimulation then
            printfn "Simulation mode, checking lglk..."
            let lglk = logLikelihood predictorsT observationsT v_p0 v_p2 v_kp d_p0 d_p2 d_kp v_bg d_bg D Sigma
            printfn "lglk %g" lglk;
            assert(abs(lglk-reference_lglk)<1e-10)
            printfn "lglk matches accpected value. Simulating..."
        else
            printfn "lglk check disabled. Simulating..."

        let testWind = getWind predictorsT v_p0 v_p2 v_kp d_p0 d_p2 d_kp

        let testWind = expandWind testWind

        let time_start = timeToTick 1500.0        
        let time_stop = timeToTick 2000.0        
        let time_by = 100
        let space_start = 0
        let space_end = 1587
    
        //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
        let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
        let dsCsv = Microsoft.Research.Science.Data.DataSet.Open("msds:csv?openMode=create&file=simulation.csv&appendMetadata=false&saveHeader=true")
        ds.IsAutocommitEnabled <- false
        dsCsv.IsAutocommitEnabled <- false
        //let windVar = ds.AddVariable<float>("windDens",[|"x";"t"|])
        let windSpVar = ds.AddVariable<float>("avgWindSpeed",[|"x";"t"|])
        //let windSpMaxVar = ds.AddVariable<float>("maxWindSpeed",[|"x";"t"|])
        let timeAxis = ds.AddVariable<int>("t",[|"t"|])
        let spaceAxis = ds.AddVariable<int>("x",[|"x"|])        

        let timeObsAxis = dsCsv.AddVariable<int>("obs_time",[|"t_obs"|])
        let obsAxis = dsCsv.AddVariable<float>("obs",[|"t_obs"|])
        let predAxis = dsCsv.AddVariable<float>("pred",[|"t_obs"|])
        let predMeanAxis = dsCsv.AddVariable<float>("predMean",[|"t_obs"|])

        printfn "Simulating for observation points..."
        let dummy =
            Table.MapToColumn "pred" ["ts";"velocity_mean"]  (fun (t:float) (v:float) ->
                let burstsAtEarth = locationState testWind (timeToTick t) (int(round(D)))
                let windAtEart = windAvg burstsAtEarth
                let v_avg_Earth,d_avg_Earth = windAtEart
                let mu = (v_bg*d_bg+v_avg_Earth*d_avg_Earth)/(d_bg+d_avg_Earth)
                        
                let n = Statistics.Normal(mu,Sigma)
                let pred = Statistics.draw rng n
                timeObsAxis.Append([|timeToTick t|]);
                obsAxis.Append([|v|]);
                predAxis.Append([|toKmPerS pred|]);
                predMeanAxis.Append([|toKmPerS mu|]);
                //printfn "%f %f %f %f" t v pred mu
                pred
                ) observationsT

        printfn "Simulated for %d time moments" (dummy.RowsCount)
        dsCsv.Commit();
        printfn "Simulating and sampling in full scale..."
        
        let simulation = simulate testWind time_start time_stop time_by space_start space_end

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
        
    *)
    0 // return an integer exit code
