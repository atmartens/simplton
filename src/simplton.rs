/***
SIMP'LTON: simulation of polysome translation
    Copyright (C) 2025 Andrew T. Martens

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License, version 3,
    as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

e-mail: simplton@logical.bio
***/

/* simplton.rs: the core logic for handling polysomes. */

use anyhow::{Result, bail};

use std::collections::VecDeque;

use rand::Rng;
use rand_distr::{StandardUniform};

use arrow2::array::*;
use arrow2::chunk::Chunk;

use std::sync::Arc;

#[inline(always)]
fn get_next_event_time(
    events : &Vec<(f32, Vec<usize>, Vec<f32>)>
) -> f32 {
    events.last()
	.map(|&(time, _, _)| time)
	.unwrap_or(f32::INFINITY)
}

/*
NOTE can also use built-in functionality for sampling distributions, but this
code is more mathematically explicit, and also maybe more performant.
 */
#[inline(always)]
fn sample_exp<R: Rng + ?Sized>(rate: f32, rng: &mut R) -> f32 {
    if rate < 0.0 {
	return f32::NAN;
    }
    else if rate == 0.0 {
	return f32::INFINITY;
    }

    let u: f32 = rng.sample(StandardUniform);
    -u.ln() / rate
}

/// Returns `Some((idx, time))` for the unblocked element
/// with the smallest sampled Exp(rate) time, or `None` if
/// everything was blocked.
fn find_next_event<R: Rng + ?Sized>(
    positions: &VecDeque<usize>,
    translation_rates: &Vec<f32>,
    blocked: &VecDeque<bool>,
    rng : &mut R
) -> Option<(usize, f32)> {
    let mut best_time = f32::INFINITY;
    let mut best_idx: Option<usize> = None;

    for (i, (&position, &is_blocked)) in
	positions
	.iter()
        .zip(blocked.iter())
        .enumerate() {
            if is_blocked {
		continue;
            }

	    // Read in the new translation rate for the next codon position
	    // +1 because P-site is pos, but A-site is rate.
	    let rate = translation_rates[position.checked_add(1).unwrap()];
	    let time = sample_exp(rate, rng);
	    
            if time < best_time {
		best_time = time;
		best_idx = Some(i);
            }
	}

    // One-liner: map an Option to a Some(index, time) pair, or a None
    best_idx.map(|i| (i, best_time))
}

pub fn run_simplton(
    start_time: f32,
    end_time: f32,
    translation_rates : &mut Vec<f32>,
    mRNA_number: u16,
    mRNA_decay_rate : f32,
    events : &mut Vec<(f32, Vec<usize>, Vec<f32>)>, // time, [pos], [rate]
    footprint_length : usize, // 10 for E. coli
    a_pos_length: usize // 4 for E. coli
) -> Result<Chunk<Arc<dyn Array>>> {
    if start_time > end_time {
	bail!("simulation cannot start after simulation end");
    }

    // Define some constants
    let delta_length : usize = footprint_length - a_pos_length;
    let max_codons = translation_rates.len() - 1;

    // RNG
    let mut rng = rand::rng();

    /* Polysome: each ribosome on the polysome is a row  */
    // Set max. capacity to num codons / footprint_length + 1; a good guess.
    let max_capacity   : usize = (max_codons / footprint_length) + 1;
    let mut identities : VecDeque<u16>   = VecDeque::with_capacity(max_capacity);
    let mut positions  : VecDeque<usize> = VecDeque::with_capacity(max_capacity);
    let mut blocked    : VecDeque<bool>  = VecDeque::with_capacity(max_capacity);
    
    /* History: each step taken is recorded */
    let mut history_identities : Vec<u16> = Vec::new();
    let mut history_positions  : Vec<u16> = Vec::new();
    let mut history_timestamps : Vec<f32> = Vec::new();

    /* Special case: initiation rate; alternates between 0.0 (when blocked)
      & true rate (when free) */
    let mut init_blocked : bool = false;

    // Has the mRNA decayed?
    let mut mRNA_decayed : bool = false;

    // Begin with an initiation event
    let mut simulation_time : f32 = start_time;

    // Keep track of largest ever ribosome id
    let mut max_id : u16 = 0;

    /*
    NOTE:
    Extensive use of checked_add() & checked_sub(), lest there be integer
    overflow or underflow. It makes the code uglier, but also more trustworthy!
    We can guarantee that any array indexing errors due to integer over/underflow
    will cause a panic.
    */
    while simulation_time < end_time {
	// 4 possible things can happen:
	// 1. A ribosome moves
	// Find earliest time & corresponding index amongst polysome
	let result = find_next_event(
	    &positions, &translation_rates, &blocked, &mut rng
	);

	// 2. A ribosome initiates
	// Inititation?
	let init_time : f32;
	
	// Initiation is blocked or mRNA has decayed
	if (init_blocked == true) || (mRNA_decayed == true) {
	    init_time = f32::INFINITY;
	}
	// Not blocked
	else {
	    init_time = sample_exp(translation_rates[0], &mut rng);
	}

	// 3. Pre-programmed (scheduled) events (changes to translation rates)
	let next_scheduled_event_time = get_next_event_time(&events);
	    
	// 4. The mRNA decays. We model this as the stop codon being inactivated
	// But don't try if already decayed!
	let mRNA_decay_time : f32;
	if mRNA_decayed == false {
	    mRNA_decay_time = sample_exp(mRNA_decay_rate, &mut rng);
	}
	else {
	    mRNA_decay_time = f32::INFINITY;
	}
	
	// Which of the 4 things actually happens?
	
	// 1. Elongation or termination for ribosome at position "index"
	// Doesn't happen if init_time happens first
	if result.is_some() &&
	    (result.unwrap().1 < init_time) &&
	    (result.unwrap().1 + simulation_time < next_scheduled_event_time) &&
	    (result.unwrap().1 < mRNA_decay_time) {
	    let result = result.unwrap();
	    let (index, dt) = (result.0 as usize, result.1);

	    // End the simulation if the delta pushes past the end time
	    if simulation_time + dt >= end_time {
		// NOTE if we ever want to read out the simulation time then it
		// makes sense to update it here, otherwise it's not necessary.
		break;
	    }

	    // Update simulation time
	    simulation_time += dt;

	    // Update historical data: timestamps
	    history_timestamps.push(simulation_time);
	    
	    // Advance the ribosome by 1 codon position
	    positions[index] = positions[index] + 1;
	    let pos = positions[index];

	    // Make sure that the ribosome didn't advance beyond the stop codon
	    assert!(pos <= max_codons);

	    // Make sure nothing weird happened: that a ribosome didn't go
	    // beyond the ribosome in front of it.
	    if index != (positions.len() - 1) {
		assert!(pos < positions[index.checked_add(1).unwrap()]);
	    }

	    // Update historical data: identities & positions
	    history_identities.push(identities[index]);
	    history_positions.push(pos as u16);	    

	    // Update any ribosome blockages
	    // But don't do it if it's already at the stop codon!
	    if pos < max_codons - 1 {
		// If it's the furthest ribosome then it can't be blocked
		if index == positions.len() - 1 {
		    // do nothing! was already unblocked
		}
		// Check if ribosome is blocked
		// NOTE: < check (as opposed to ==) is crucial to handle edge
		// cases around a ribosome which recently initiated!
		else if pos + footprint_length
		    < positions[index.checked_add(1).unwrap()] {
		    blocked[index] = false;
		}
		// Not the only ribosome, and failed the collision check.
		else {
		    blocked[index] = true;
		}
	    }
	    
	    // Unblock translation initiation?
	    // NOTE: special A site positioning calculation for initiation.
	    if index == 0 {
		if pos > delta_length {
		    init_blocked = false;
		}
		// Special case: if a termination ribosome unblocks initiation
		// Can only happen on very short mRNAs which only fit 1 ribosome
		// at a time.
		else if (pos == max_codons) && (pos <= delta_length) {
		    init_blocked = false;
		}
	    }
	    // Unblock previous ribosome?
	    else {
		if pos == positions[index.checked_sub(1).unwrap()]
		    + footprint_length + 1 {
			blocked[index.checked_sub(1).unwrap()] = false;
		    }
	    }
		
	    // Termination?
	    if pos == max_codons {
		// Remove the ribosome from the polysome
		identities.pop_back().unwrap();
		positions.pop_back().unwrap();
		blocked.pop_back().unwrap();
	    }
	    // Does the ribosome which just advanced become blocked?
	    else if (index != positions.len() - 1)
		&& (pos + footprint_length
		    == positions[index.checked_add(1).unwrap()]) {
		    blocked[index] = true;
	    }
	}
	// 2. Initiation took place instead of elongation or scheduled event AND
	// mRNA has not decayed
	else if (init_time + simulation_time < next_scheduled_event_time) &&
	    (init_time < mRNA_decay_time) {
	    // End the simulation if the delta pushes past the end time
	    if simulation_time + init_time >= end_time {
		break;
	    }
	    
	    // Update simulation time
	    simulation_time += init_time;

	    // Update historical data: timestamps
	    history_timestamps.push(simulation_time);

	    let id = max_id;
	    max_id = max_id.checked_add(1).unwrap();

	    // Update historical data: identities & positions
	    history_identities.push(id as u16);
	    history_positions.push(0);

	    identities.push_front(id);
	    positions.push_front(0);

	    // Single ribosome
	    if positions.len() == 1 {
		blocked.push_front(false);
	    }
	    // polysome, but not blocked
	    else if positions[1] > delta_length {
		blocked.push_front(false);
	    }
	    // Blocked
	    else {
		blocked.push_front(true);
	    }
	    
	    // And now set the initiation codon to blocked
	    init_blocked = true;
	}
	// 3. If we get this far then it must be a pre-programmed event (rate
	// change)
	else if next_scheduled_event_time < simulation_time + mRNA_decay_time {
	    // Special case: pre-programmed event happens next, but it actually
	    // is beyond the simulation end time, so do nothing & end the
	    // simulation
	    if next_scheduled_event_time >= end_time {
		break;
	    }
	    else {
		let (now, codons, tr_rates) = events.pop().unwrap();

		// Update translation_rates
		for (p, r) in codons.iter().zip(tr_rates.iter()) {
		    translation_rates[*p] = *r;
		}
		
		// Advance the simulation
		simulation_time = now;
	    }
	}
	// Special case: all times are infinity, and therefore equal to one
	// another. The simulation needs to terminate.
	else if result.is_none() &&
	    (init_time == f32::INFINITY) &&
	    (next_scheduled_event_time == f32::INFINITY) &&
	    (mRNA_decay_time ==  f32::INFINITY) {
	    break;
	}
	// 4. Last possibility: mRNA decay
	else {
	    if simulation_time + mRNA_decay_time >= end_time {
		break;
	    }
	    mRNA_decayed = true;
	    simulation_time += mRNA_decay_time;
	}
    }

    /* Simulation Finished */
    
    // Return a "chunk" of arrays
    let num_columns = 4;
    let mut arrays: Vec<Arc<dyn Array>> = Vec::with_capacity(num_columns);

    let history_length = history_identities.len();

    // Include the mRNA simulation number as a column
    let mRNA_number_list: Vec<u16> = vec![mRNA_number; history_length];
    let mRNA_number_array = PrimitiveArray::from_vec(mRNA_number_list);
    arrays.push(mRNA_number_array.arced());

    // Now insert each of the vectors into "arrays"
    arrays.push(PrimitiveArray::from_vec(history_identities).arced());
    arrays.push(PrimitiveArray::from_vec(history_positions).arced());
    arrays.push(PrimitiveArray::from_vec(history_timestamps).arced());

    let trajectories = Chunk::try_new(arrays)?;

    Ok(trajectories)
}

