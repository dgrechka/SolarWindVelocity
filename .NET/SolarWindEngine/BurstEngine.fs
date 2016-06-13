module BurstEngine

(**
#   Solar Wind Computational model

## Overview

The model is descrete. It presents time as descrete ticks

*)

type Time = int // world ticks


(**
The space of the world is descrete offset from the sun
*)

type Location = int // descrete steps from the Sun

let descretize (l:float) : Location =
    int(System.Math.Round(l))

(**

The world consists of a number of Sun wind bursts emmited from the Sun
The velocity in the model is number of steps per model tick.
*)

type Burst = {
    EmergenceTime: Time    
    Velocity: float //portion of max velocity: between 0.0 and 1.0
    Density: float
    }

type Wind = Burst list

(**
At each time the world state is charactirized by the position of all bursts
*)

type WorldState = Map<Location,Set<Burst>> // At each location can be any number of bursts

let getBurstPosition burst time =
    let pos = burst.Velocity * float(time-burst.EmergenceTime)
    descretize(System.Math.Round(pos))

let worldState wind time =
    let rec fillState wind (state:WorldState) =
        match wind with
        |   burst::tail  ->
            let pos = getBurstPosition burst time
            let burstsAtPos =
                if Map.containsKey pos state then
                    state.[pos]
                else
                    Set.empty
            let burstAtPosAppended = Set.add burst burstsAtPos           
            fillState tail (Map.add pos burstAtPosAppended state)
        |   []  ->
            state
    fillState wind Map.empty

let windAvg (bursts:Set<Burst>) =
        let folder acc burst =
            let acc_v,acc_d = acc
            let d = burst.Density
            acc_v+(burst.Velocity*d),acc_d+d
        let accumulated = Set.fold folder (0.0,0.0) bursts
        let acc_v,acc_d = accumulated
        acc_v/acc_d,acc_d

let windAtDistance (world_state:WorldState) d =
    if world_state.ContainsKey d then
        windAvg world_state.[d]
    else
        0.0,0.0
        

let sampleVelocity wind time max_x =    
    let state = worldState wind time        
    Array.init max_x (windAtDistance state) |> Array.map fst


//outer list is time, inner list is space
let simulate wind start_t stop_t left_x right_x =
    let N = stop_t - start_t + 1    
    Array.Parallel.init N ( fun i ->
        printf "."
        let t = start_t+i        
        sampleVelocity wind t (right_x+1) |> Array.skip (left_x-1)
        ) |> List.ofArray