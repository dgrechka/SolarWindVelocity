﻿module BurstEngine

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
            if pos < 0 then //optimizing. does not account bursts that are not yet emerged
                fillState tail state
            else
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

let locationState wind time location =
    let rec fillState wind (state:Set<Burst>) =        
        match wind with
        |   burst::tail  ->
            let pos = getBurstPosition burst time
            if pos = location then //looking for particular location
                fillState tail (Set.add burst state) //mixing in current burst                
            else                                
                fillState tail state //not interesting location, skipping
        |   []  ->
            state
    fillState wind Set.empty

let expandWind wind = //interpolates the bursts linearly, so burst appear every tick
    let interpolate b1 b2 = //returns reversed list
        let t1 = b1.EmergenceTime
        let t2 = b2.EmergenceTime
        let steps = t2-t1
        assert(steps >= 1)        
        let v_step = (b2.Velocity-b1.Velocity)/float(steps)
        let d_step = (b2.Density-b1.Density)/float(steps)
        [
            for i=1 to steps do
                yield {
                    EmergenceTime=t2-i
                    Velocity=b2.Velocity - v_step*float(i)
                    Density=b2.Density - d_step*float(i)
                }
        ]
    let rec expand expanded unexpanded =
        match unexpanded with
        |   h1::h2::tail ->
            expand (List.append (interpolate h1 h2) expanded) (h2::tail)
        |   h1::[] -> h1::expanded
        |   []  -> expanded
    expand [] wind |> List.rev

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
    Array.init (max_x+1) (windAtDistance state) |> Array.map fst


//outer list is time, inner list is space
let simulate wind start_t stop_t left_x right_x =
    let N = stop_t - start_t + 1    
    Array.Parallel.init N ( fun i ->
        printf "."
        let t = start_t+i        
        sampleVelocity wind t right_x |> Array.skip left_x
        ) |> List.ofArray