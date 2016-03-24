module app

open PulseEngine

open Angara
open Angara.Data
open Angara.Statistics
open Angara.Filzbach

open MathNet.Numerics.LinearAlgebra

let getWind (predictorsT: Table) V_b S_0 a b c P_w= 
    let generatePulse t area = 
        if area < S_0 then
            None
        else        
            let v = a*area*area+b*area+c                
            assert(a*S_0*S_0+b*S_0+c = V_b)
            let pulse = {
                Power=1.0;
                EmergenceTime=t;
                Velocity=v*60.0*60.0;
                A=v-V_b;
                Kernel=Kernels.GaussianExt P_w
            }
            Some(pulse)
    Table.Map ["t";"CHarea"] generatePulse predictorsT |> Seq.choose (fun a -> a) |> List.ofSeq

let mutable bestLglk = System.Double.NegativeInfinity

let get2ndPolyCoefs x0 x_offset1 x_offset2 y0 y_offset1 y_offset2 =
    let x1 = x0+x_offset1
    let x2 = x1+x_offset2    
    let m = matrix [[ x0*x0; x0; 1.0]
                    [ x1*x1; x1; 1.0]
                    [ x2*x2; x2; 1.0]]
    let y1=y0+y_offset1
    let y2=y1+y_offset2
    let r_part = vector [y0; y1; y2]
    let solution = m.Solve r_part
    assert(Vector.length solution = 3)
    (solution.[0],solution.[1],solution.[2])

let logLikelihood (predictorsT: Table) (observationsT: Table) (p: Parameters) =                
    let V_b = p.GetValue "V_b" //"background" velocity
    let V_1 = p.GetValue "V_1" //first node's v component
    let V_2 = p.GetValue "V_2" //second node's v component
    let S_0 = p.GetValue "S_0" //minimal CH area that enables slope
    let S_1 = p.GetValue "S_1" //first node's s component
    let S_2 = p.GetValue "S_2" //second node's s component

    //Distance to Earth
    //varies between 0.9832898912 and 1.0167103335 AU. 
    //which is 1.47098074 to 1.52097701 ×10^8 km
    let D = p.GetValue "D" //distance to Earth
    let Sigma = p.GetValue "Sigma" //observation noise
    let P_w = p.GetValue "P_w" //Pulse width. around pulse distance traveled in half an hour?

    let (a,b,c) = get2ndPolyCoefs S_0 S_1 S_2 V_b V_1 V_2

    let wind = getWind predictorsT V_b S_0 a b c P_w
    
    let data = Table.Map ["t";"velocity"] (fun (t:float) (v:float) -> (t,v)) observationsT

    let current_lglk p =
        let t,v = p
        let maxVelocity = WindAMaxForTime wind t
        let predAtEarth = maxVelocity D
        let mu = V_b + predAtEarth //base level plus pulses        
        if System.Double.IsNaN mu then
            log improbable
        else
            let distribution = Statistics.Normal(mu,Sigma)
            let lglk = log_pdf distribution v
            //printfn "mu:%g\tv:%g\tsigma:%g\tlglk:%g" mu v Sigma lglk
            lglk    
    let res = data |> Seq.filter (fun t -> not(System.Double.IsNaN(snd t))) |> Array.ofSeq |> Array.map current_lglk |> Seq.sum         
    if res>bestLglk then
        printfn "\nlglk:%g\tV: %g\t%g\t%g\tS: %g\t%g\t%g\tP_w:%g\tSig:%g\tD:%g" res V_b (V_b+V_1) (V_b+V_1+V_2) S_0 (S_0+S_1) (S_0+S_1+S_2) P_w Sigma D
        bestLglk <- res
    else
        printf "."
    res
        

[<EntryPoint>]
let main argv =     
    //tests
    let (a,b,c) = get2ndPolyCoefs -1.0 1.0 1.0 2.0 -1.0 1.0
    assert(a=1.0)
    assert(b=0.0)
    assert(c=1.0)

    let (a,b,c) = get2ndPolyCoefs 0.0 1.0 1.0 2.0 -1.0 5.0
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

    let ch_data = Table.Load @"..\..\..\..\1.ChAreaApproximated.csv"
    let obs = Table.Load @"..\..\..\..\SampleData\ace_swepam_2015_1h.csv"
    let obs = Table.MapToColumn "t" ["hours since 01.01.2015  1:00:00"]  (fun t -> t+1.0) obs
    
    let mutable counter = 0
    let obs =
        if thinnObs then
            Table.Filter ["t"]  (fun (dummy:float) ->            
            counter <- counter+1
            counter%5=0            
            ) obs
        else obs

    let predictorsT = Table.Filter ["t"] (fun t -> t > 40.0 ) ch_data
    let observationsT = Table.Filter ["t"] (fun t -> t > 150.0 && t < 2140.0) obs

    printfn "Predictor values count %d" predictorsT.RowsCount
    printfn "Observation values count %d" observationsT.RowsCount    
    
    let eps = System.Double.Epsilon

    let estimate predictorsT observationsT =
        let parameters = Parameters.Empty
                            .Add("V_b",[|350.0|],eps,500.0,isLog=true,delay=0)
                            .Add("V_1",[|37.6|],eps,200.0,isLog=false,delay=0)
                            .Add("V_2",[|37.6|],eps,200.0,isLog=false,delay=0)
                            .Add("S_0",[|2.1|],eps,10.0,isLog=true,delay=0)
                            .Add("S_1",[|1.0|],eps,10.0,isLog=false,delay=0)
                            .Add("S_2",[|1.0|],eps,10.0,isLog=false,delay=0)
                            .Add("D",[|1.83843e+08|],1.40e+8,2.00e+8,isLog=false,delay=0)
                            .Add("Sigma",[| 71.4027|],eps, 500.0,isLog=true,delay=0)
                            .Add("P_w",[|3.77509e+06|],300.0*60.0*10.0, 300.0*60.0*210.0,isLog=false,delay=0)        
        Sampler.runmcmc(parameters, logLikelihood predictorsT observationsT, 5000, 1000,rng=rng)
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
        let V_b,V_1,V_2,S_0,S_1,S_2,P_w,Sigma,D = 350.0,40.06,27.8,0.18,2.78,1.04,3.76868e+06,70.7075,1.78556e+08

        let (a,b,c) = get2ndPolyCoefs S_0 S_1 S_2 V_b V_1 V_2

        let testWind = getWind predictorsT V_b S_0 a b c P_w

        let time_start = 0.0
        let time_step = 0.2
        let time_stop = 500.0
        let space_step = 0.5e+6
        let space_start = 0.0 
        let space_end = 1.5e+8

        let simulation = simulate testWind time_start time_stop time_step space_start space_end space_step

        //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
        let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
        ds.IsAutocommitEnabled <- false
        let windVar = ds.AddVariable<float>("windDens",[|"x";"t"|])
        let windSpVar = ds.AddVariable<float>("avgWindSpeed",[|"x";"t"|])
        let windSpMaxVar = ds.AddVariable<float>("maxWindSpeed",[|"x";"t"|])
        let timeAxis = ds.AddVariable<float>("t",[|"t"|])
        let spaceAxis = ds.AddVariable<float>("x",[|"x"|])        

        let timeObsAxis = ds.AddVariable<float>("obs_time",[|"t_obs"|])
        let obsAxis = ds.AddVariable<float>("obs",[|"t_obs"|])
        let predAxis = ds.AddVariable<float>("pred",[|"t_obs"|])
        let predMeanAxis = ds.AddVariable<float>("predMean",[|"t_obs"|])

        let dummy =
            Table.MapToColumn "pred" ["t";"velocity"]  (fun (t:float) (v:float) ->
                let mu = V_b + WindAvgAForTime testWind t D
                let n = Statistics.Normal(mu,Sigma)
                let pred = Statistics.draw rng n
                timeObsAxis.Append([|t|]);
                obsAxis.Append([|v|]);
                predAxis.Append([|pred|]);
                predMeanAxis.Append([|mu|]);
                printfn "%f %f %f %f" t v pred mu
                pred
                ) observationsT

        printfn "simulated for %d time moments" (dummy.RowsCount)

        List.iteri (fun i sim_step_data ->
            let sample_vals,sample_speed,sample_max_speed = sim_step_data
            let sample_speed_km_s = List.map (fun a -> a + V_b) sample_speed
            let t = time_start+float(i)*time_step                            
            windVar.Append(List.toArray sample_vals,"t")
            windSpVar.Append(List.toArray sample_speed_km_s,"t")
            windSpMaxVar.Append(List.toArray sample_max_speed,"t")
            timeAxis.Append([|t|]);
            )
            simulation
        let spaceAxisData = Array.init (int((space_end-space_start)/space_step)) (fun i -> space_start+float(i)*space_step)
        spaceAxis.Append(spaceAxisData)

        ds.Commit()        
    
    0 // return an integer exit code
