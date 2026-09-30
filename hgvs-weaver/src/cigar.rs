use crate::error::HgvsError;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq)]
pub struct CigarOp {
    pub op: char,
    pub len: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cigar {
    pub ops: Vec<CigarOp>,
}

impl Cigar {
    pub fn ref_len(&self) -> i32 {
        self.ops
            .iter()
            .filter(|op| "=MXDN".contains(op.op))
            .map(|op| op.len)
            .sum()
    }

    pub fn tgt_len(&self) -> i32 {
        self.ops
            .iter()
            .filter(|op| "=MXI".contains(op.op))
            .map(|op| op.len)
            .sum()
    }
}

impl FromStr for Cigar {
    type Err = HgvsError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ops = Vec::new();
        let mut current_num = String::new();
        for c in s.chars() {
            if c.is_ascii_digit() {
                current_num.push(c);
            } else {
                let len = if current_num.is_empty() {
                    1
                } else {
                    current_num.parse::<i32>().map_err(|_| {
                        HgvsError::Other(format!("Invalid number in CIGAR string: {}", current_num))
                    })?
                };
                current_num.clear();
                ops.push(CigarOp { op: c, len });
            }
        }
        Ok(Cigar { ops })
    }
}

impl std::fmt::Display for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for op in &self.ops {
            if op.len != 1 {
                write!(f, "{}", op.len)?;
            }
            write!(f, "{}", op.op)?;
        }
        Ok(())
    }
}

/// Where a position on one side of an alignment falls on the other side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Projected {
    /// The position is aligned to a base: its counterpart's position.
    Base(i32),
    /// The position has no counterpart (an `I` seen from the transcript, a
    /// `D` seen from the genome): it lies between positions `at - 1` and `at`
    /// on the other side. `run` is the whole run of such positions it belongs
    /// to, half-open, on its own side.
    Gap { at: i32, run: (i32, i32) },
    /// The position is inside an `N` (an intron): the nearer flanking base
    /// and the signed distance to it, as an intronic offset.
    Intronic { base: i32, offset: i32 },
}

#[derive(Debug, Clone)]
pub struct CigarMapper {
    pub cigar: Cigar,
    pub ref_pos: Vec<i32>,
    pub tgt_pos: Vec<i32>,
}

impl CigarMapper {
    pub fn new(cigar_str: &str) -> Result<Self, HgvsError> {
        let cigar = Cigar::from_str(cigar_str)?;
        let mut ref_pos = Vec::with_capacity(cigar.ops.len() + 1);
        let mut tgt_pos = Vec::with_capacity(cigar.ops.len() + 1);

        let mut ref_cur = 0;
        let mut tgt_cur = 0;

        for op in &cigar.ops {
            ref_pos.push(ref_cur);
            tgt_pos.push(tgt_cur);

            if "=MXDN".contains(op.op) {
                ref_cur += op.len;
            }
            if "=MXI".contains(op.op) {
                tgt_cur += op.len;
            }
        }
        ref_pos.push(ref_cur);
        tgt_pos.push(tgt_cur);

        Ok(CigarMapper {
            cigar,
            ref_pos,
            tgt_pos,
        })
    }

    pub fn ref_len(&self) -> i32 {
        *self.ref_pos.last().unwrap()
    }

    pub fn tgt_len(&self) -> i32 {
        *self.tgt_pos.last().unwrap()
    }

    /// Where reference position `pos` falls on the target. A `D` base is a
    /// `Gap`; an `N` base is `Intronic`.
    pub fn map_ref_to_tgt(&self, pos: i32, strict_bounds: bool) -> Result<Projected, HgvsError> {
        self.map_internal(&self.ref_pos, &self.tgt_pos, pos, strict_bounds)
    }

    /// Where target position `pos` falls on the reference. An `I` base is a
    /// `Gap`.
    pub fn map_tgt_to_ref(&self, pos: i32, strict_bounds: bool) -> Result<Projected, HgvsError> {
        self.map_internal(&self.tgt_pos, &self.ref_pos, pos, strict_bounds)
    }

    fn map_internal(
        &self,
        from_pos: &[i32],
        to_pos: &[i32],
        pos: i32,
        strict_bounds: bool,
    ) -> Result<Projected, HgvsError> {
        let last_pos = *from_pos.last().unwrap();
        if strict_bounds && (pos < 0 || pos > last_pos) {
            return Err(HgvsError::Other(
                "Position is beyond the bounds of sequence".to_string(),
            ));
        }

        let mut pos_i = 0;
        for i in 0..self.cigar.ops.len() {
            if pos < from_pos[i + 1] {
                pos_i = i;
                break;
            }
            pos_i = i;
        }

        match self.cigar.ops[pos_i].op {
            '=' | 'M' | 'X' => Ok(Projected::Base(to_pos[pos_i] + (pos - from_pos[pos_i]))),
            'D' | 'I' => Ok(Projected::Gap {
                at: to_pos[pos_i],
                run: (from_pos[pos_i], from_pos[pos_i + 1]),
            }),
            'N' => {
                if pos - from_pos[pos_i] < from_pos[pos_i + 1] - pos {
                    Ok(Projected::Intronic {
                        base: to_pos[pos_i] - 1,
                        offset: pos - from_pos[pos_i] + 1,
                    })
                } else {
                    Ok(Projected::Intronic {
                        base: to_pos[pos_i],
                        offset: -(from_pos[pos_i + 1] - pos),
                    })
                }
            }
            op => Err(HgvsError::Other(format!("Unsupported CIGAR op: {}", op))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Projected::{Base, Gap, Intronic};

    #[test]
    fn test_cigarmapper() {
        let cigar_str = "3=2N=X=3N=I=D=";
        let cm = CigarMapper::new(cigar_str).unwrap();

        assert_eq!(cm.ref_len(), 15);
        assert_eq!(cm.tgt_len(), 10);

        // ref to tgt
        assert_eq!(cm.map_ref_to_tgt(0, true).unwrap(), Base(0));
        assert_eq!(
            cm.map_ref_to_tgt(3, true).unwrap(),
            Intronic { base: 2, offset: 1 }
        );
        assert_eq!(
            cm.map_ref_to_tgt(4, true).unwrap(),
            Intronic {
                base: 3,
                offset: -1
            }
        );
        assert_eq!(cm.map_ref_to_tgt(5, true).unwrap(), Base(3));
        assert_eq!(cm.map_ref_to_tgt(6, true).unwrap(), Base(4));
        // ref 12 lands in op 8 (=) because op 7 (I) doesn't consume ref
        assert_eq!(cm.map_ref_to_tgt(12, true).unwrap(), Base(8));
        // ref 13 is the D: the target lacks it, and it lies between tgt 8 and 9
        assert_eq!(
            cm.map_ref_to_tgt(13, true).unwrap(),
            Gap {
                at: 9,
                run: (13, 14)
            }
        );
        assert_eq!(cm.map_ref_to_tgt(14, true).unwrap(), Base(9));

        // tgt to ref
        assert_eq!(cm.map_tgt_to_ref(0, true).unwrap(), Base(0));
        assert_eq!(cm.map_tgt_to_ref(3, true).unwrap(), Base(5));
        assert_eq!(cm.map_tgt_to_ref(4, true).unwrap(), Base(6));
        // tgt 7 is the I: the reference lacks it, and it lies between ref 11 and 12
        assert_eq!(
            cm.map_tgt_to_ref(7, true).unwrap(),
            Gap {
                at: 12,
                run: (7, 8)
            }
        );
        // tgt 8 lands in op 8 (=)
        assert_eq!(cm.map_tgt_to_ref(8, true).unwrap(), Base(12));
    }

    #[test]
    fn test_cigarmapper_strict_bounds() {
        let cigar_str = "3=2N=X=3N=I=D=";
        let cm = CigarMapper::new(cigar_str).unwrap();

        assert!(cm.map_ref_to_tgt(-1, true).is_err());
        assert!(cm.map_ref_to_tgt(16, true).is_err());

        assert_eq!(cm.map_ref_to_tgt(0, true).unwrap(), Base(0));
        assert_eq!(cm.map_ref_to_tgt(-1, false).unwrap(), Base(-1));
        assert_eq!(cm.map_ref_to_tgt(15, true).unwrap(), Base(10));
        assert_eq!(cm.map_ref_to_tgt(14, false).unwrap(), Base(9));
    }

    #[test]
    fn a_gap_at_the_edge_of_an_alignment_is_placed_at_its_end() {
        // Soft-clipped starts and ends come as leading and trailing I ops.
        let cm = CigarMapper::new("5I20=").unwrap();
        assert_eq!(
            cm.map_tgt_to_ref(0, true).unwrap(),
            Gap { at: 0, run: (0, 5) }
        );
        assert_eq!(
            cm.map_tgt_to_ref(4, true).unwrap(),
            Gap { at: 0, run: (0, 5) }
        );
        assert_eq!(cm.map_tgt_to_ref(5, true).unwrap(), Base(0));
        let cm = CigarMapper::new("20=5I").unwrap();
        assert_eq!(cm.map_tgt_to_ref(19, true).unwrap(), Base(19));
        assert_eq!(
            cm.map_tgt_to_ref(20, true).unwrap(),
            Gap {
                at: 20,
                run: (20, 25)
            }
        );
        assert_eq!(
            cm.map_tgt_to_ref(24, true).unwrap(),
            Gap {
                at: 20,
                run: (20, 25)
            }
        );
    }

    #[test]
    fn test_cigar_parsing() {
        let cigar_str = "3=2N=X=3N=I=D=";
        let cigar = Cigar::from_str(cigar_str).unwrap();
        assert_eq!(cigar.ops.len(), 11);
        assert_eq!(cigar.ops[0], CigarOp { op: '=', len: 3 });
        assert_eq!(cigar.to_string(), cigar_str);
    }
}
