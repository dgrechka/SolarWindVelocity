module PulseEngine

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


let WindAvgVelocityForTime wind time = 
    let SheftedKernelFun pulse time = 
        let delta = (time-pulse.EmergenceTime)*pulse.Velocity
        let f t = (pulse.Kernel (t-delta))
        f
    let avg t =
        let portions = List.map (fun pulse -> SheftedKernelFun pulse time t) wind
        let total = List.sum portions
        let total_back = 1.0/total
        let portions_velocities = List.zip portions (List.map (fun pulse -> pulse.Velocity) wind)
        let velocity = List.map (fun pair -> let p,v = pair in p*v*p*total_back) portions_velocities |> List.sum
        velocity
    avg

let sample (f:float -> float) left right (by:float) =
    let N = int((right-left)/by)
    List.init N (fun i -> f (left+by*float(i)))

//outer list is time, inner list is space
let simulate wind start_t stop_t (step_t:float) left_x right_x by_x =
    //Checking kernels
    List.iter (fun pulse -> assert(pulse.Kernel(0.0)=1.0)) wind

    let N = int((stop_t-start_t)/step_t)
    List.init N ( fun i ->
        let t = start_t+step_t*float(i)
        sample (WindFuncForTime wind t) left_x right_x by_x,
        sample (WindAvgVelocityForTime wind t) left_x right_x by_x
        )

type Kernels () =
    static member Step t = if t < -1.0 || t >= 1.0 then 0.0 else 1.0
    static member Gaussian t = exp(-t*t/2.0)
