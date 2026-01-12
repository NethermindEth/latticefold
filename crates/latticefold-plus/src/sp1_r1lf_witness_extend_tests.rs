#![cfg(all(test, feature = "we_gate"))]

use crate::sp1_r1lf_witness_extend::orig_num_vars_from_counts;

#[test]
fn test_orig_num_vars_from_counts() {
    assert_eq!(orig_num_vars_from_counts(10, 3).unwrap(), 7);
    assert!(orig_num_vars_from_counts(3, 10).is_err());
}

