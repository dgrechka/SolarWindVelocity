type Time = float

type KernelFunc = float -> float //real function that has a peek at 0.0

type Pulse = {
    EmergenceTime: Time
    Power: double //1.0 is default
    Velocity: float
    Kernel: KernelFunc
    }

type Wind = Pulse list

let PulseFuncForTime pulse time =
    let delta = (time-pulse.EmergenceTime)*pulse.Velocity
    let f t = (pulse.Kernel (t-delta))*pulse.Power
    f

let WindFuncForTime wind time =
    let sum t =
        let folder (acc:float) pulse =
            acc+ (PulseFuncForTime pulse time t)
        List.fold folder 0.0 wind
    sum

let sample (f:float -> float) left right (by:float) =
    let N = int((right-left)/by)
    List.init N (fun i -> f (left+by*float(i)))

//outer list is time, inner list is space
let simulate wind start_t stop_t (step_t:float) left_x right_x by_x =
    let N = int((stop_t-start_t)/step_t)
    List.init N ( fun i ->
        let t = start_t+step_t*float(i)
        sample (WindFuncForTime wind t) left_x right_x by_x
        )

[<EntryPoint>]
let main argv = 
    let stepKernel t = if t < -1.0 || t >= 1.0 then 0.0 else 1.0
    let gaussianKernel t = exp(-t*t/2.0)
    
    let pulse1 = {
        EmergenceTime= 0.0;
        Power= 1.0
        Velocity= 1.0;
        Kernel= gaussianKernel
    }

    let pulse2 = { //slow
        EmergenceTime= 10.0;
        Power= 0.5
        Velocity= 0.5;
        Kernel= gaussianKernel
    }

    let pulse3 = { //fast
        EmergenceTime= 20.0;
        Power= 2.0
        Velocity= 2.0;
        Kernel= gaussianKernel
    }
    
    let testWind = [pulse1; pulse2; pulse3]
    
    let time_start = 0.0
    let time_step = 0.1
    let space_step = 0.1
    let space_start = 0.0
    let space_end = 90.0
        

    let simulation = simulate testWind time_start 60.0 time_step space_start space_end space_step

    //Dumping the data to NetCDF using http://research.microsoft.com/en-us/downloads/ccf905f6-34c6-4845-892e-a5715a508fa3/
    let ds = Microsoft.Research.Science.Data.DataSet.Open("msds:nc?openMode=create&file=simulation.nc")
    ds.IsAutocommitEnabled <- false
    let windVar = ds.AddVariable<float>("wind",[|"x";"t"|])
    let timeAxis = ds.AddVariable<float>("t",[|"t"|])
    let spaceAxis = ds.AddVariable<float>("x",[|"x"|])    

    List.iteri (fun i sample_vals ->
        let t = time_start+float(i)*time_step
        windVar.Append(List.toArray sample_vals,"t")
        timeAxis.Append([|t|]);
        )
        simulation
    let spaceAxisData = Array.init (int((space_end-space_start)/space_step)) (fun i -> space_start+float(i)*space_step)
    spaceAxis.Append(spaceAxisData)

    ds.Commit()        
    0 // return an integer exit code
