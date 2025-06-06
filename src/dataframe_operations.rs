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

/* dataframe_operations.rs: calculate statistics on trajectories using Polars 
dataframes. */

use polars::prelude::*;
use polars::datatypes::DataType::{UInt16, Float32, Int64};
use polars::series::ops::NullBehavior::Ignore;

/// NOTE: Warmup time: useful if only want statistics after this time
/// otherwise default to 0.0 to include all activity.

// NOTE uses a dataframe with complete step times
pub fn calculate_mean_step_times(
    df_step_times: &DataFrame
) -> PolarsResult<DataFrame> {
    let sort_opts = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true);

    let result = df_step_times
	.clone()
	.lazy()
	.group_by(["pos"])
	.agg([
	    col("step time").mean()
	])
	.sort(["pos"], sort_opts)
	.collect()?;
    
    Ok(result)
}

// FIXME getting a null pos with null step time
// Where does it come from? Just remove it?
pub fn calculate_complete_step_times(
    df_trajectories: &DataFrame,
    df_mrna_start_times : &DataFrame,
    warmup_time: f32,
) -> PolarsResult<DataFrame> {
    // Makes the code easier to read to pre-define these
    let sort_opts_double = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending_multi([false, false]);

    let sort_opts_triple = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending_multi([false, false, false]);

    let union_args = UnionArgs {
	parallel : true,
	rechunk : true,
	to_supertypes: false,
	diagonal : false,
	from_partitioned_ds : false,
	maintain_order : true
    };

    
    let step_times = df_trajectories
	.clone()
	.lazy()
	.filter(col("timestamp").gt(warmup_time))
	.sort(["mRNA", "ribosome", "pos"], sort_opts_triple.clone())
        .with_column(
	    col("timestamp")
		.diff(1, Ignore)
		.over([col("mRNA"), col("ribosome")])
		.alias("step time")
	)
	.drop(["timestamp"])
	.sort(["mRNA", "ribosome"], sort_opts_double.clone());
    
    let init_times = df_trajectories
	.clone()
	.lazy()
        .filter(col("pos").eq(0))
	.filter(col("timestamp").gt(warmup_time))
	.sort(["ribosome", "timestamp"], sort_opts_double.clone())
        .with_column(
	    col("timestamp")
		.diff(1, Ignore)
		.over([col("mRNA")])
		.alias("step time")
	)
	.drop(["pos", "timestamp"])	      
        .sort(["mRNA"], Default::default());

    let init_zero = df_trajectories
	.clone()
	.lazy()
        .filter(
	    col("ribosome")
		.eq(0)
		.and(col("pos")
		     .eq(0))
	)
        .sort(["mRNA"], Default::default());

    let df_mrna_start_times = df_mrna_start_times
	.clone()
	.lazy()
	.sort(["mRNA"], Default::default())
	.drop(["mRNA"]);

    let concat_zero_and_start = concat_lf_horizontal(
	&[init_zero, df_mrna_start_times],
	union_args.clone()
    )?
        .with_column(
	    (col("timestamp") - col("start time"))
		.alias("step time")
	)
	      .drop(["timestamp", "start time"]);

    let init_times_no_pioneer = init_times
        .filter(col("ribosome").neq(0))
        .with_column(
	    lit(0_u16)
	    	.alias("pos")
		.cast(DataType::UInt16)
	    // NOTE/FIXME why do we need to cast??
	)
        .select(&[col("mRNA"), col("ribosome"), col("pos"), col("step time")]);

    
    let init_times_yes_pioneer = concat(
	&[concat_zero_and_start, init_times_no_pioneer],
	union_args
    )?
	.sort(["mRNA", "ribosome"], sort_opts_double.clone());

    let step_times_no_start = step_times
        .filter(col("pos").neq(0));

    let complete_step_times = concat(
	&[init_times_yes_pioneer, step_times_no_start],
	union_args
    )?
	.sort(["mRNA", "ribosome"], sort_opts_double.clone())
	// FIXME Getting a null pos + step time? For now just drop that row.
	.drop_nulls(None)
	.collect()?;

    Ok(complete_step_times)
}

pub fn calculate_final_occupancy(
    df_trajectories: &DataFrame,
    num_codons: u16,
    num_mRNA : u16,
    footprint_length: usize,
    a_site_length: usize
) -> PolarsResult<DataFrame> {
    let footprint_length = footprint_length as u16;
    let a_site_length = a_site_length as u16;
    let delta_length : u16 = footprint_length - a_site_length;

    let sort_opts = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending(false);
    
    let result = df_trajectories
	.clone()
	.lazy()
	.drop(["timestamp"])
	.filter(
	    // Must not have yet terminated
	    col("pos")
		.max()
		.over(["mRNA", "ribosome"])
		.lt(num_codons - 1)
	)
    // Offset footprint
	.with_column(col("pos") + lit(delta_length))
	.group_by(["mRNA", "ribosome"])
	.agg(
	    [col("pos").top_k(lit(footprint_length)).alias("pos")]
	)
	.explode(["pos"])
    // Don't count occupancy beyond the stop codon
	.filter(col("pos").lt_eq(num_codons))
	.group_by(["pos"])
	.agg([col("pos").count().cast(Float32).alias("count")])
	.with_column(
	    (col("count") / lit(num_mRNA as f32)).alias("fractional occupancy")
	)
	.drop(["count"])
	.sort(["pos"], sort_opts)
	.collect()?;
    
    Ok(result)
}

// TODO: arbitrary time?
pub fn calculate_final_nascent_peptide_length_distribution(
    df_trajectories : &DataFrame,
    num_codons : u16
) -> PolarsResult<DataFrame> {
    let sort_opts = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending(false);

    let result = df_trajectories
	.clone()
	.lazy()
	.filter(
	// Must not have yet terminated
	col("pos")
	    .max()
	    .over(["mRNA", "ribosome"])
	    .lt(num_codons - 1)
    )
	.drop(["timestamp", "mRNA", "ribosome"])
	.group_by(["pos"])
	.agg([col("pos").count().alias("count")])
	.sort(["pos"], sort_opts)
	.collect()?;
    
    Ok(result)
}

pub fn calculate_protein_production_by_mrna(
    df_trajectories: &DataFrame,
    num_codons: u16,
    warmup_time: f32
) -> PolarsResult<DataFrame> {
    let sort_opts_double = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending_multi([false, false]);
    
    let protein_production = df_trajectories
	.clone()
	.lazy()
        .filter(col("pos").eq(num_codons-1))
	.filter(col("timestamp").gt(warmup_time))
	.sort(["mRNA", "timestamp"], sort_opts_double.clone())
        .with_column(lit(1).alias("unit"))
        .select(&[
            col("mRNA"),
            col("timestamp"),
            col("unit")
		.cum_sum(false)
		.over([col("mRNA")])
		.alias("total"),
        ])
	.sort(["mRNA", "timestamp"], sort_opts_double)
        .collect()?;

    Ok(protein_production)
}

pub fn calculate_total_protein_production(
    df_trajectories: &DataFrame,
    num_codons: u16,
    warmup_time: f32
) -> PolarsResult<DataFrame> {
    let protein_production = df_trajectories
	.clone()
	.lazy()
        .filter(col("pos").eq(num_codons - 1))
	.filter(col("timestamp").gt(warmup_time))
	.sort(["timestamp"], Default::default())
	.with_column(lit(1_i64).alias("unit"))
	.with_column(
	    col("unit")
		.cum_sum(false)
		.alias("total")
	)
	.select(&[col("mRNA"),col("timestamp"),	col("total")])
	.sort(["timestamp"], Default::default())
	.collect()?;
    
    Ok(protein_production)
}

pub fn calculate_pulse_chase(
    df_trajectories: DataFrame,
    pulse_time: f32,
    chase_time: f32,
    label_positions: Vec<u16>,
    num_codons: u16
) -> PolarsResult<DataFrame> {
    let sort_opts_triple = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending_multi([false, false, false]);

    let l_p = Series::new("label_positions".into(), label_positions.clone());

    let labeled = df_trajectories
	.clone()
	.lazy()
	.filter(
	    col("pos")
		.max()
		.over(&["mRNA", "ribosome"])
		.eq(num_codons - 1)
	).filter(
	    col("pos")
		.is_in(lit(l_p))
	).filter(
	    col("timestamp")
		.gt(pulse_time)
	).with_column(
	    lit(1_u16)
		.alias("radioactivity")
	).with_column(
	    col("timestamp")
		.max()
		.over(&["mRNA", "ribosome"])
		.alias("exit time")
	).filter(
	    col("timestamp")
		.lt(chase_time)
	).sort(
	    ["mRNA", "ribosome", "pos"],
	    sort_opts_triple
	).with_column(
	    col("radioactivity")
		.cum_sum(false)
		.over(&["mRNA", "ribosome"])
		.alias("radioactivity")
	).sort(
	    ["exit time"], Default::default()
	).with_column(
	    col("radioactivity")
		.cum_sum(false)
		.cast(Float32)
		.alias("radioactivity")
	).with_column(
	    (
		col("radioactivity") /
		    col("radioactivity").max() *
		    lit(label_positions.len() as f32)
	    ).alias("radioactivity")
	).with_column(
	    (col("exit time") - col("exit time").min())
		.alias("relative time")
	).collect()?;

	Ok(labeled)
}

pub fn calculate_synthesis_times(
    df_trajectories: &DataFrame,
    num_codons: u16,
    warmup_time: f32
) -> PolarsResult<DataFrame> {
    let sort_opts_double = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_nulls_last(true)
	.with_order_descending_multi([false, false]);

    let df_differences = df_trajectories
	.clone()
	.lazy()
    // Must have finished translation
	.filter(
	    col("pos")
		.max()
		.over(["mRNA", "ribosome"])
		.eq(num_codons - 1)
	)
        .group_by(&[col("mRNA"), col("ribosome")])
        .agg([
            col("timestamp").min().alias("init time"),
            col("timestamp").max().alias("term time"),
        ])
	.filter(col("init time").gt(warmup_time))
	.filter(col("term time").gt(warmup_time))
        .with_column(
	    (col("term time") - col("init time")).alias("synthesis time")
	)
        .drop(["init time", "term time"])
    // NOTE Arguably instead want to sort by radioactivity, exit time, or
    // relative time so that plotting doesn't require sorting later?
	.sort(["mRNA", "ribosome"], sort_opts_double)
        .collect()?;

    Ok(df_differences)
}

pub fn calculate_final_average_ribosome_spacing(
    df_trajectories: &DataFrame,
    num_codons: u16
) -> PolarsResult<DataFrame> {
    let final_ribosomes = df_trajectories
	.clone()
	.lazy()
	.drop(["timestamp"])
	.with_column(
	    col("pos")
		.max()
		.over(["mRNA", "ribosome"])
		.alias("max_pos")
	)
    // Must not have yet terminated
	.filter(col("max_pos").lt(num_codons - 1))
    // Must be the furthest it goes
	.filter(col("pos").eq(col("max_pos")))
	.drop(["max_pos"])
	.with_column((col("ribosome") + lit(1)).alias("next_ribosome"));

    let spacings = final_ribosomes.clone().join(
	final_ribosomes.clone(),
	&[ col("mRNA"), col("ribosome")],
	&[ col("mRNA"), col("next_ribosome")],
	JoinArgs::new(JoinType::Inner).with_suffix(Some("_next".into()))
    )
	.with_column((col("pos_next") - col("pos")).alias("spacing"))
	.drop(["pos"])
	.group_by(["mRNA"])
	.agg([col("spacing").mean().alias("average spacing")])
	.collect()?;

    Ok(spacings)
}

/*
How to calculate collisions

Take trajectories, groupby 1. mRNA # 2. ribosome # 3. position
    -> each such unique combination yields a unique timestamp.

Then do a self-join twice:
First, on (ribosome -1, pos-1) (if it exists: exception is either initiation, or last ribosome to initiate)
Second, on (ribosome, pos+1) (if it exists: exception is either termination, or end of simulation)
And run 1 comparison:
is ribosome-1,pos-1 timestamp < or >= ribosome,pos+1 timestamp?

If ribosome-1,pos-1 < ribosome,pos+1 then it's a collision. Else it's not.

Weird situation with equality: a bit ambigous as whether to call it a collision,
but I'd argue that hypothetically it would /not/ be, because if the forward ribosome
advances in unison with the one before, then by definition it did not block it!

On the other hand, if two ribosomes advance in unison, that would kind of
suggest that the forwards ribosome /is/ rate-limiting... ?

In practice, it's infintessimally unlikely that two events would happen at the same time
in f32 space, and in practice I don't think the random number generator can actually
generate 0.0 time events, so the act of picking a single event in the case of a tie
actually makes the scenario impossible in practice. So it's more a matter of convention
than handling an edge case.
 */

pub fn calculate_collision_frequency(
    df_trajectories: &DataFrame,
    num_codons : u16,
    footprint_length : usize,
    warmup_time: f32
) -> PolarsResult<DataFrame> {
    let footprint_length = footprint_length as u16;
    
    let sort_opts = SortMultipleOptions::default()
	.with_maintain_order(false)
	.with_multithreaded(true)
	.with_nulls_last(true)
	.with_order_descending(false);

    let result = df_trajectories
	.clone()
	.lazy()
	.filter(col("timestamp").gt(warmup_time))
	.with_column(col("pos").cast(Int64))
	.with_column((col("pos") + lit(1)).alias("pos_next"));

    let result = result
	.clone()
	.join(
	    result.clone(),
	    &[ col("mRNA"), col("ribosome"), col("pos_next") ],
	    &[ col("mRNA"), col("ribosome"), col("pos") ],
	    JoinArgs::new(JoinType::Left).with_suffix(Some("_advance".into())))
	.drop(["pos_next", "pos_next_advance"])
	.with_columns([
	    (col("ribosome") + lit(1)).alias("rib_neighbor"),
	    (col("pos") - lit(footprint_length)).alias("pos_neighbor")
	]);

    let result = result.lazy();

    let result = result
	.clone()
	.join(
	    result.clone(),
	    &[ col("mRNA"), col("rib_neighbor"), col("pos_neighbor") ],
	    &[ col("mRNA"), col("ribosome"), col("pos") ],
	    JoinArgs::new(JoinType::Left).with_suffix(Some("_collide".into())))
	.drop(["rib_neighbor", "pos_neighbor", "rib_neighbor_collide",
	       "pos_neighbor_collide"])
	.with_columns([
	    col("timestamp_advance").fill_null(f32::INFINITY),
	    col("timestamp_collide").fill_null(f32::INFINITY)
	])
	.with_column(
	    col("timestamp_collide")
		.lt(col("timestamp_advance"))
		.cast(UInt16)
		.alias("is_collision")
	)
    // Remove pos == stop codon since by definition a ribosome can't engage
    // in collisions after it's been removed from the mRNA.
	.filter(col("pos").lt(num_codons - 1))
	.group_by(["pos"])
	.agg([
	    col("is_collision").sum().cast(Float32).alias("collision_count"),
	    col("is_collision").count().cast(Float32).alias("total_events"),
	])
	.with_column(
	    ((col("collision_count") * lit(100.0_f32)) /
	     col("total_events")).alias("collision_percentage")
	)
	.sort(["pos"], sort_opts)
	.collect()?;

    Ok(result)
}

