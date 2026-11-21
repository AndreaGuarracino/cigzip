mod agc_index;
mod faidx;
mod gdb_index;
mod sequence_index;

// Re-export public API
pub use sequence_index::{collect_sequence_paths, SequenceIndex};

/// Reused per-thread buffers larger than this are freed after use.
pub const BUF_KEEP: usize = 4 << 20;
