//! Suffix array (libsais) + LCP (Kasai) over the integer mirror encoding.
//!
//! libsais builds the SA over `sa_input`. LCP is computed by Kasai over the
//! original (untouched) `sa_input` rather than chaining libsais' PLCP, because
//! the large-alphabet path takes the text by `&mut` and may clobber it. Unique
//! sentinels make every suffix distinct and cap each LCP at the first special,
//! reproducing gt's special-aware LCP. lcp[0] = 0; lcp[r] = lcp(SA[r-1], SA[r]).

use libsais::{suffix_array::AlphabetSize, SuffixArrayConstruction};

/// Returns (suffix_array, lcp) over all `n = sa_input.len()` positions.
pub fn sa_lcp(sa_input: &[i32], k: i32) -> (Vec<i32>, Vec<i32>) {
    let mut text = sa_input.to_vec(); // libsais large-alphabet path mutates text
    let mut c = SuffixArrayConstruction::for_text_mut(&mut text)
        .in_owned_buffer32()
        .single_threaded();
    // SAFETY: k is max(sa_input) + 1 by construction (dense 0..k-1).
    unsafe {
        c = c.with_alphabet_size(AlphabetSize::new(k));
    }
    let res = c.run().expect("libsais suffix array");
    let sa: Vec<i32> = res.suffix_array().to_vec();
    drop(res); // release the &mut text borrow
    let lcp = kasai(sa_input, &sa);
    (sa, lcp)
}

/// Kasai's O(n) LCP. lcp[rank] = lcp(SA[rank-1], SA[rank]); lcp[0] = 0.
fn kasai(text: &[i32], sa: &[i32]) -> Vec<i32> {
    let n = sa.len();
    let mut rank = vec![0i32; n];
    for (i, &s) in sa.iter().enumerate() {
        rank[s as usize] = i as i32;
    }
    let mut lcp = vec![0i32; n];
    let mut h = 0usize;
    for i in 0..n {
        if rank[i] > 0 {
            let j = sa[(rank[i] - 1) as usize] as usize;
            while i + h < n && j + h < n && text[i + h] == text[j + h] {
                h += 1;
            }
            lcp[rank[i] as usize] = h as i32;
            if h > 0 {
                h -= 1;
            }
        } else {
            h = 0;
        }
    }
    lcp
}

#[cfg(test)]
mod tests {
    use super::*;

    // Cross-check libsais SA + Kasai LCP against a brute-force suffix sort on a
    // small integer text with unique sentinels.
    #[test]
    fn sa_lcp_matches_naive() {
        let text: Vec<i32> = vec![0, 1, 2, 3, 0, 1, 2, 3, 4]; // ACGTACGT + sentinel
        let (sa, lcp) = sa_lcp(&text, 5);
        let n = text.len();
        let mut naive: Vec<usize> = (0..n).collect();
        naive.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
        assert_eq!(sa.iter().map(|&x| x as usize).collect::<Vec<_>>(), naive);
        // verify lcp against direct computation
        for r in 1..n {
            let (a, b) = (sa[r - 1] as usize, sa[r] as usize);
            let mut l = 0;
            while a + l < n && b + l < n && text[a + l] == text[b + l] {
                l += 1;
            }
            assert_eq!(lcp[r] as usize, l, "lcp mismatch at rank {r}");
        }
        assert_eq!(lcp[0], 0);
    }
}
