mod faidx;
mod sequence_index;

#[cfg(feature = "agc")]
mod agc_index;

// Re-export public API
pub use sequence_index::{collect_sequence_paths, SequenceIndex};
