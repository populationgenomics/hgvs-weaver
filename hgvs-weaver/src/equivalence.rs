use crate::analogous_edit::{project_aa_variant, project_na_variant, SparseReference};
use crate::data::{DataProvider, IdentifierKind, TranscriptSearch};
use crate::error::HgvsError;
use crate::mapper::VariantMapper;
use crate::structs::{GVariant, IntervalSpdi, NaEdit, PVariant, SequenceVariant, Variant};
use crate::utils::decompose_aa;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquivalenceLevel {
    /// Identical notation after basic normalization.
    Identity,
    /// Biologically identical but different notation (e.g., ins vs dup).
    Analogous,
    /// Definitively different edits/outcomes.
    Different,
    /// Missing data or unsupported variant type for comparison.
    Unknown,
}

impl EquivalenceLevel {
    pub fn is_equivalent(&self) -> bool {
        matches!(self, Self::Identity | Self::Analogous)
    }
}

// Migrated to analogous_edit.rs

fn strand_aware_edit(edit: &NaEdit, strand: crate::data::Strand) -> NaEdit {
    if strand == crate::data::Strand::Minus {
        edit.reverse_complement()
    } else {
        edit.clone()
    }
}

pub struct VariantEquivalence<'a> {
    pub hdp: &'a dyn DataProvider,
    pub searcher: &'a dyn TranscriptSearch,
    pub mapper: VariantMapper<'a>,
}

impl<'a> VariantEquivalence<'a> {
    pub fn new(hdp: &'a dyn DataProvider, searcher: &'a dyn TranscriptSearch) -> Self {
        VariantEquivalence {
            hdp,
            searcher,
            mapper: VariantMapper::new(hdp),
        }
    }

    pub fn equivalent(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<bool, HgvsError> {
        Ok(self.equivalent_level(var1, var2)?.is_equivalent())
    }

    pub fn equivalent_level(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        // Expand gene symbols if present
        let vars1 = self.expand_if_gene_symbol(var1)?;
        let vars2 = self.expand_if_gene_symbol(var2)?;

        for v1 in &vars1 {
            for v2 in &vars2 {
                let lvl = self.equivalent_level_single(v1, v2)?;
                if lvl.is_equivalent() {
                    return Ok(lvl);
                }
            }
        }
        Ok(EquivalenceLevel::Different)
    }

    fn equivalent_level_single(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<EquivalenceLevel, HgvsError> {
        // 1. Strict Check (after normalization)
        if self.normalize_format(&var1.to_string()) == self.normalize_format(&var2.to_string()) {
            return Ok(EquivalenceLevel::Identity);
        }

        // 2. Build and Merge Sparse References
        let s1 = self.get_ref_for_variant(var1);
        let s2 = self.get_ref_for_variant(var2);
        let mut merged = s1;
        if let Err(_) = merged.merge(&s2) {
            return Ok(EquivalenceLevel::Different); // Inconsistent references
        }

        // 3. Project and Compare Outcomes
        match (var1, var2) {
            (SequenceVariant::Protein(p1), SequenceVariant::Protein(p2)) => {
                let pos1_opt = &p1.posedit.pos;
                let pos2_opt = &p2.posedit.pos;

                // Handle global identity (p.=) vs positional variant
                let (start1, end1, edit1, start2, end2, edit2) = match (pos1_opt, pos2_opt) {
                    (Some(pos1), Some(pos2)) => {
                        let s1 = pos1.start.base.to_index().0;
                        let e1 = self.get_effective_end(p1, s1);
                        let s2 = pos2.start.base.to_index().0;
                        let e2 = self.get_effective_end(p2, s2);
                        (s1, e1, &p1.posedit.edit, s2, e2, &p2.posedit.edit)
                    }
                    (Some(pos1), None) if p2.posedit.edit.is_identity() => {
                        let s1 = pos1.start.base.to_index().0;
                        let e1 = self.get_effective_end(p1, s1);
                        // Synthesize p2 interval to match p1
                        (s1, e1, &p1.posedit.edit, s1, e1, &p2.posedit.edit)
                    }
                    (None, Some(pos2)) if p1.posedit.edit.is_identity() => {
                        let s2 = pos2.start.base.to_index().0;
                        let e2 = self.get_effective_end(p2, s2);
                        // Synthesize p1 interval to match p2
                        (s2, e2, &p1.posedit.edit, s2, e2, &p2.posedit.edit)
                    }
                    _ => {
                        // Fallback to cross-type comparison logic which handles non-projected cases
                        if self.are_equivalent_single(var1, var2)? {
                            return Ok(EquivalenceLevel::Analogous);
                        }
                        return Ok(EquivalenceLevel::Different);
                    }
                };

                let min_pos = start1.min(start2);
                let max_pos = end1.max(end2);

                let res1 = project_aa_variant(edit1, start1, end1, min_pos, max_pos, &merged)
                    .trim_at_stop();
                let res2 = project_aa_variant(edit2, start2, end2, min_pos, max_pos, &merged)
                    .trim_at_stop();

                let is_analogous = res1.is_analogous_to(&res2);

                if is_analogous {
                    return Ok(EquivalenceLevel::Analogous);
                }
            }
            (SequenceVariant::Coding(c1), SequenceVariant::Coding(c2)) => {
                if let (Some(pos1), Some(pos2)) = (&c1.posedit.pos, &c2.posedit.pos) {
                    let mut i1 = pos1.spdi_interval(&c1.ac, self.hdp)?;
                    let mut i2 = pos2.spdi_interval(&c2.ac, self.hdp)?;

                    let t1 = self.hdp.get_transcript(&c1.ac, None)?;
                    let edit1 = strand_aware_edit(&c1.posedit.edit, t1.strand);

                    let t2 = self.hdp.get_transcript(&c2.ac, None)?;
                    let edit2 = strand_aware_edit(&c2.posedit.edit, t2.strand);

                    if matches!(c1.posedit.edit, NaEdit::Ins { .. }) && pos1.end.is_some() {
                        // An insertion between two flanking bases is anchored at the
                        // lower genomic index for projection.
                        let (p, _, ac) = i1;
                        i1 = (p, p + 1, ac);
                    }
                    if matches!(c2.posedit.edit, NaEdit::Ins { .. }) && pos2.end.is_some() {
                        // An insertion between two flanking bases is anchored at the
                        // lower genomic index for projection.
                        let (p, _, ac) = i2;
                        i2 = (p, p + 1, ac);
                    }

                    let (start1, end1, _) = i1;
                    let (start2, end2, _) = i2;

                    let min_pos = start1.min(start2).saturating_sub(2);
                    let max_pos = end1.max(end2) + 2;

                    let res1 =
                        project_na_variant(&edit1, start1, end1 - 1, min_pos, max_pos - 1, &merged);
                    let res2 =
                        project_na_variant(&edit2, start2, end2 - 1, min_pos, max_pos - 1, &merged);

                    if res1.is_analogous_to(&res2) {
                        return Ok(EquivalenceLevel::Analogous);
                    }
                }
            }
            (SequenceVariant::NonCoding(n1), SequenceVariant::NonCoding(n2)) => {
                if let (Some(pos1), Some(pos2)) = (&n1.posedit.pos, &n2.posedit.pos) {
                    let mut i1 = pos1.spdi_interval(&n1.ac, self.hdp)?;
                    let mut i2 = pos2.spdi_interval(&n2.ac, self.hdp)?;

                    let t1 = self.hdp.get_transcript(&n1.ac, None)?;
                    let edit1 = strand_aware_edit(&n1.posedit.edit, t1.strand);

                    let t2 = self.hdp.get_transcript(&n2.ac, None)?;
                    let edit2 = strand_aware_edit(&n2.posedit.edit, t2.strand);

                    if matches!(n1.posedit.edit, NaEdit::Ins { .. }) && pos1.end.is_some() {
                        // An insertion between two flanking bases is anchored at the
                        // lower genomic index for projection.
                        let (p, _, ac) = i1;
                        i1 = (p, p + 1, ac);
                    }
                    if matches!(n2.posedit.edit, NaEdit::Ins { .. }) && pos2.end.is_some() {
                        // An insertion between two flanking bases is anchored at the
                        // lower genomic index for projection.
                        let (p, _, ac) = i2;
                        i2 = (p, p + 1, ac);
                    }

                    let (start1, end1, _) = i1;
                    let (start2, end2, _) = i2;

                    let min_pos = start1.min(start2).saturating_sub(2);
                    let max_pos = end1.max(end2) + 2;

                    let res1 =
                        project_na_variant(&edit1, start1, end1 - 1, min_pos, max_pos - 1, &merged);
                    let res2 =
                        project_na_variant(&edit2, start2, end2 - 1, min_pos, max_pos - 1, &merged);

                    if res1.is_analogous_to(&res2) {
                        return Ok(EquivalenceLevel::Analogous);
                    }
                }
            }
            _ => {
                // Fallback to existing logic for cross-type comparison
                if self.are_equivalent_single(var1, var2)? {
                    if self.is_cross_type_identity(var1, var2) {
                        return Ok(EquivalenceLevel::Identity);
                    }
                    return Ok(EquivalenceLevel::Analogous);
                }
            }
        }

        Ok(EquivalenceLevel::Different)
    }

    fn is_cross_type_identity(&self, var1: &SequenceVariant, var2: &SequenceVariant) -> bool {
        match (var1, var2) {
            (SequenceVariant::Coding(vc), SequenceVariant::Protein(vp))
            | (SequenceVariant::Protein(vp), SequenceVariant::Coding(vc)) => {
                if let Ok(vp_generated) = self.mapper.c_to_p(vc, Some(&vp.ac)) {
                    vp_generated.to_string() == vp.to_string()
                } else {
                    false
                }
            }
            (SequenceVariant::Genomic(vg), SequenceVariant::Coding(vc))
            | (SequenceVariant::Coding(vc), SequenceVariant::Genomic(vg)) => {
                if let Ok(tx) = self.hdp.get_transcript(&vc.ac, None) {
                    if let Ok(vg_generated) = self
                        .mapper
                        .c_to_g(vc, Some(tx.reference_accession.as_str()))
                    {
                        vg_generated.to_string() == vg.to_string()
                    } else {
                        false
                    }
                } else {
                    false
                }
            }
            (SequenceVariant::Genomic(vg), SequenceVariant::NonCoding(vn))
            | (SequenceVariant::NonCoding(vn), SequenceVariant::Genomic(vg)) => {
                if let Ok(tx) = self.hdp.get_transcript(&vn.ac, None) {
                    if let Ok(vg_generated) = self
                        .mapper
                        .n_to_g(vn, Some(tx.reference_accession.as_str()))
                    {
                        vg_generated.to_string() == vg.to_string()
                    } else {
                        false
                    }
                } else {
                    false
                }
            }
            (SequenceVariant::NonCoding(vn), SequenceVariant::Protein(vp))
            | (SequenceVariant::Protein(vp), SequenceVariant::NonCoding(vn)) => {
                if let Ok(tx) = self.hdp.get_transcript(&vn.ac, None) {
                    if let Ok(vg_generated) = self
                        .mapper
                        .n_to_g(vn, Some(tx.reference_accession.as_str()))
                    {
                        if let Ok(c_variants) = self.mapper.g_to_c_all(&vg_generated, self.searcher)
                        {
                            for vc in c_variants {
                                if let Ok(vp_generated) = self.mapper.c_to_p(&vc, Some(&vp.ac)) {
                                    if vp_generated.to_string() == vp.to_string() {
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                }
                false
            }
            _ => false,
        }
    }

    fn get_effective_end(&self, vp: &PVariant, start: i32) -> i32 {
        let mut end = vp.posedit.pos.as_ref().map_or(start, |pos| {
            pos.end.as_ref().map_or(start, |e| e.base.to_index().0)
        });

        if let crate::structs::AaEdit::Repeat { ref_: Some(s), .. } = &vp.posedit.edit {
            let len = s.len() as i32;
            if end - start + 1 < len {
                end = start + len - 1;
            }
        }
        end
    }

    fn get_ref_for_variant(&self, var: &SequenceVariant) -> SparseReference {
        let mut s = SparseReference::new();
        match var {
            SequenceVariant::Protein(vp) => {
                if let Ok(seq) = self
                    .mapper
                    .refs
                    .reference(&vp.ac, crate::data::IdentifierType::ProteinAccession)
                    .whole()
                {
                    if let Ok(aas) = decompose_aa(&seq) {
                        for (i, aa) in aas.iter().enumerate() {
                            let _ = s.set(i as i32, aa.to_string());
                        }
                    }
                }
            }
            SequenceVariant::Coding(vc) => {
                if let Some(pos) = &vc.posedit.pos {
                    if let Ok((start, end, spdi_ac)) = pos.spdi_interval(&vc.ac, self.hdp) {
                        if let (Ok(s0), Ok(e0)) = (usize::try_from(start), usize::try_from(end)) {
                            if let Ok(seq) = self
                                .mapper
                                .refs
                                .reference(&spdi_ac, crate::data::IdentifierType::GenomicAccession)
                                .slice(s0, e0)
                            {
                                let _ = s.set(start, seq);
                            }
                        }
                    }
                }
            }
            _ => {}
        }
        s
    }

    fn expand_if_gene_symbol(
        &self,
        var: &SequenceVariant,
    ) -> Result<Vec<SequenceVariant>, HgvsError> {
        let ac = var.ac();

        // Use DataProvider to determine if this is a symbol or an accession.
        let id_type = self.hdp.get_identifier_type(ac)?;

        if id_type == crate::data::IdentifierType::GeneSymbol {
            let target_kind = match var {
                SequenceVariant::Protein(_) => IdentifierKind::Protein,
                SequenceVariant::Coding(_)
                | SequenceVariant::NonCoding(_)
                | SequenceVariant::Rna(_) => IdentifierKind::Transcript,
                _ => IdentifierKind::Genomic,
            };

            // Try symbol expansion.
            let accessions =
                self.hdp
                    .get_symbol_accessions(ac, IdentifierKind::Genomic, target_kind)?;

            if !accessions.is_empty() {
                let mut expanded = Vec::new();
                for (ac_type, new_ac) in accessions {
                    // Only expand into compatible types.
                    let is_compatible = match (var, ac_type) {
                        (
                            SequenceVariant::Protein(_),
                            crate::data::IdentifierType::ProteinAccession,
                        ) => true,
                        (
                            SequenceVariant::Coding(_)
                            | SequenceVariant::NonCoding(_)
                            | SequenceVariant::Rna(_),
                            crate::data::IdentifierType::TranscriptAccession,
                        ) => true,
                        (
                            SequenceVariant::Genomic(_) | SequenceVariant::Mitochondrial(_),
                            crate::data::IdentifierType::GenomicAccession,
                        ) => true,
                        // Allow g. on transcripts if specifically provided
                        (
                            SequenceVariant::Genomic(_),
                            crate::data::IdentifierType::TranscriptAccession,
                        ) => true,
                        _ => false,
                    };

                    if is_compatible {
                        let mut v = var.clone();
                        match &mut v {
                            SequenceVariant::Genomic(v_g) => v_g.ac = new_ac,
                            SequenceVariant::Coding(v_c) => v_c.ac = new_ac,
                            SequenceVariant::Protein(v_p) => v_p.ac = new_ac,
                            SequenceVariant::Mitochondrial(v_m) => v_m.ac = new_ac,
                            SequenceVariant::NonCoding(v_n) => v_n.ac = new_ac,
                            SequenceVariant::Rna(v_r) => v_r.ac = new_ac,
                        }
                        expanded.push(v);
                    }
                }
                if !expanded.is_empty() {
                    return Ok(expanded);
                }
            }
        }

        // Default: return as-is.
        Ok(vec![var.clone()])
    }

    fn are_equivalent_single(
        &self,
        var1: &SequenceVariant,
        var2: &SequenceVariant,
    ) -> Result<bool, HgvsError> {
        match (var1, var2) {
            // Nucleotide vs Nucleotide
            (SequenceVariant::Genomic(v1), SequenceVariant::Genomic(v2)) => {
                self.n_vs_n_equivalent(v1, v2)
            }
            (SequenceVariant::Coding(v1), SequenceVariant::Coding(v2)) => {
                self.n_vs_n_equivalent_c(v1, v2)
            }
            (SequenceVariant::NonCoding(v1), SequenceVariant::NonCoding(v2)) => {
                self.n_vs_n_equivalent_n(v1, v2)
            }

            (SequenceVariant::Genomic(v1), SequenceVariant::Coding(v2)) => {
                self.g_vs_c_equivalent(v1, v2)
            }
            (SequenceVariant::Coding(v1), SequenceVariant::Genomic(v2)) => {
                self.g_vs_c_equivalent(v2, v1)
            }

            (SequenceVariant::Genomic(v1), SequenceVariant::NonCoding(v2)) => {
                self.g_vs_n_equivalent(v1, v2)
            }
            (SequenceVariant::NonCoding(v1), SequenceVariant::Genomic(v2)) => {
                self.g_vs_n_equivalent(v2, v1)
            }

            (SequenceVariant::Coding(v1), SequenceVariant::NonCoding(v2)) => {
                self.c_vs_n_equivalent(v1, v2)
            }
            (SequenceVariant::NonCoding(v1), SequenceVariant::Coding(v2)) => {
                self.c_vs_n_equivalent(v2, v1)
            }

            // Nucleotide vs Protein
            (SequenceVariant::Genomic(v1), SequenceVariant::Protein(v2)) => {
                self.g_vs_p_equivalent(v1, v2)
            }
            (SequenceVariant::Protein(v1), SequenceVariant::Genomic(v2)) => {
                self.g_vs_p_equivalent(v2, v1)
            }
            (SequenceVariant::Coding(v1), SequenceVariant::Protein(v2)) => {
                self.c_vs_p_equivalent(v1, v2)
            }
            (SequenceVariant::Protein(v1), SequenceVariant::Coding(v2)) => {
                self.c_vs_p_equivalent(v2, v1)
            }
            (SequenceVariant::NonCoding(v1), SequenceVariant::Protein(v2)) => {
                self.n_vs_p_equivalent(v1, v2)
            }
            (SequenceVariant::Protein(v1), SequenceVariant::NonCoding(v2)) => {
                self.n_vs_p_equivalent(v2, v1)
            }

            // Protein vs Protein
            (SequenceVariant::Protein(v1), SequenceVariant::Protein(v2)) => {
                self.p_vs_p_equivalent(v1, v2)
            }

            // Fallback
            _ => {
                if var1.ac() == var2.ac() && var1.coordinate_type() == var2.coordinate_type() {
                    return Ok(self.normalize_format(&var1.to_string())
                        == self.normalize_format(&var2.to_string()));
                }
                Ok(false)
            }
        }
    }

    fn normalize_format(&self, s: &str) -> String {
        // Strip parentheses and replace '?' with 'X' (unknown amino acid, analogous to Xaa)
        let mut s = s.replace(['(', ')'], "").replace('?', "X");
        // Normalize 3-letter AA codes to 1-letter
        let replacements = [
            ("Ala", "A"),
            ("Arg", "R"),
            ("Asn", "N"),
            ("Asp", "D"),
            ("Cys", "C"),
            ("Gln", "Q"),
            ("Glu", "E"),
            ("Gly", "G"),
            ("His", "H"),
            ("Ile", "I"),
            ("Leu", "L"),
            ("Lys", "K"),
            ("Met", "M"),
            ("Phe", "F"),
            ("Pro", "P"),
            ("Ser", "S"),
            ("Thr", "T"),
            ("Trp", "W"),
            ("Tyr", "Y"),
            ("Val", "V"),
            ("Asx", "B"),
            ("Glx", "Z"),
            ("Xaa", "X"),
            ("Ter", "*"),
        ];
        for (from, to) in replacements {
            s = s.replace(from, to);
        }
        s
    }

    fn n_vs_n_equivalent(&self, v1: &GVariant, v2: &GVariant) -> Result<bool, HgvsError> {
        let nv1 = self
            .mapper
            .normalize_variant(SequenceVariant::Genomic(v1.clone()))?;
        let nv2 = self
            .mapper
            .normalize_variant(SequenceVariant::Genomic(v2.clone()))?;

        let s1 = self.normalize_format(&nv1.to_string());
        let s2 = self.normalize_format(&nv2.to_string());

        Ok(s1 == s2)
    }

    fn n_vs_n_equivalent_c(
        &self,
        v1: &crate::structs::CVariant,
        v2: &crate::structs::CVariant,
    ) -> Result<bool, HgvsError> {
        let tx1 = self.hdp.get_transcript(&v1.ac, None)?;
        let tx2 = self.hdp.get_transcript(&v2.ac, None)?;
        let g1 = self
            .mapper
            .c_to_g(v1, Some(tx1.reference_accession.as_str()))?;
        let g2 = self
            .mapper
            .c_to_g(v2, Some(tx2.reference_accession.as_str()))?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn n_vs_n_equivalent_n(
        &self,
        v1: &crate::structs::NVariant,
        v2: &crate::structs::NVariant,
    ) -> Result<bool, HgvsError> {
        let tx1 = self.hdp.get_transcript(&v1.ac, None)?;
        let tx2 = self.hdp.get_transcript(&v2.ac, None)?;
        let g1 = self
            .mapper
            .n_to_g(v1, Some(tx1.reference_accession.as_str()))?;
        let g2 = self
            .mapper
            .n_to_g(v2, Some(tx2.reference_accession.as_str()))?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn g_vs_c_equivalent(
        &self,
        vg: &crate::structs::GVariant,
        vc: &crate::structs::CVariant,
    ) -> Result<bool, HgvsError> {
        let g2 = self.mapper.c_to_g(vc, Some(&vg.ac))?;
        self.n_vs_n_equivalent(vg, &g2)
    }

    fn g_vs_n_equivalent(
        &self,
        vg: &crate::structs::GVariant,
        vn: &crate::structs::NVariant,
    ) -> Result<bool, HgvsError> {
        let g2 = self.mapper.n_to_g(vn, Some(&vg.ac))?;
        self.n_vs_n_equivalent(vg, &g2)
    }

    fn c_vs_n_equivalent(
        &self,
        vc: &crate::structs::CVariant,
        vn: &crate::structs::NVariant,
    ) -> Result<bool, HgvsError> {
        let tx = self.hdp.get_transcript(&vc.ac, None)?;
        let ref_ac = tx.reference_accession;
        let g1 = self.mapper.c_to_g(vc, Some(ref_ac.as_str()))?;
        let g2 = self.mapper.n_to_g(vn, Some(ref_ac.as_str()))?;
        self.n_vs_n_equivalent(&g1, &g2)
    }

    fn n_vs_p_equivalent(
        &self,
        vn: &crate::structs::NVariant,
        vp: &crate::structs::PVariant,
    ) -> Result<bool, HgvsError> {
        let tx = self.hdp.get_transcript(&vn.ac, None)?;
        let ref_ac = tx.reference_accession;
        let vg = self.mapper.n_to_g(vn, Some(ref_ac.as_str()))?;
        self.g_vs_p_equivalent(&vg, vp)
    }

    fn g_vs_p_equivalent(
        &self,
        vg: &crate::structs::GVariant,
        vp: &crate::structs::PVariant,
    ) -> Result<bool, HgvsError> {
        let c_variants = self.mapper.g_to_c_all(vg, self.searcher)?;
        for vc in c_variants {
            if self.c_vs_p_equivalent(&vc, vp)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn c_vs_p_equivalent(
        &self,
        vc: &crate::structs::CVariant,
        vp: &crate::structs::PVariant,
    ) -> Result<bool, HgvsError> {
        let vp_generated = self.mapper.c_to_p(vc, Some(&vp.ac))?;
        Ok(self.normalize_format(&vp_generated.to_string())
            == self.normalize_format(&vp.to_string()))
    }

    fn p_vs_p_equivalent(&self, v1: &PVariant, v2: &PVariant) -> Result<bool, HgvsError> {
        if v1.ac == v2.ac {
            if self.normalize_format(&v1.to_string()) == self.normalize_format(&v2.to_string()) {
                return Ok(true);
            }

            if let (Some(pos1), Some(pos2)) = (&v1.posedit.pos, &v2.posedit.pos) {
                let start1 = pos1.start.base.to_index().0;
                let end1 = self.get_effective_end(v1, start1);

                let start2 = pos2.start.base.to_index().0;
                let end2 = self.get_effective_end(v2, start2);

                let min_pos = start1.min(start2).saturating_sub(2);
                let max_pos = end1.max(end2) + 2;

                let mut sref = self.get_ref_for_variant(&SequenceVariant::Protein(v1.clone()));
                let sref2 = self.get_ref_for_variant(&SequenceVariant::Protein(v2.clone()));
                sref.merge(&sref2)?;

                let res1 = crate::analogous_edit::project_aa_variant(
                    &v1.posedit.edit,
                    start1,
                    end1,
                    min_pos,
                    max_pos,
                    &sref,
                );
                let res2 = crate::analogous_edit::project_aa_variant(
                    &v2.posedit.edit,
                    start2,
                    end2,
                    min_pos,
                    max_pos,
                    &sref,
                );

                return Ok(res1.is_analogous_to(&res2));
            }
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::{GenomicPos, TranscriptPos};
    use crate::data::{ExonData, IdentifierKind, IdentifierType, TranscriptData};

    struct MockDataProvider;
    impl DataProvider for MockDataProvider {
        fn get_transcript(
            &self,
            ac: &str,
            _ref_ac: Option<&str>,
        ) -> Result<TranscriptData, HgvsError> {
            if ac == "NM_000123.4" {
                Ok(TranscriptData {
                    ac: "NM_000123.4".to_string(),
                    gene: "ABC".to_string(),
                    cds_start_index: Some(TranscriptPos(0)),
                    cds_end_index: Some(TranscriptPos(19)),
                    strand: crate::data::Strand::Plus,
                    reference_accession: "NC_000001.11".to_string(),
                    exons: vec![ExonData {
                        transcript_start: TranscriptPos(0),
                        transcript_end: TranscriptPos(19),
                        reference_start: GenomicPos(0),
                        reference_end: GenomicPos(19),
                        alt_strand: crate::data::Strand::Plus,
                        cigar: "20M".to_string(),
                    }],
                })
            } else {
                Err(HgvsError::ValidationError("Not found".into()))
            }
        }
        fn get_seq(
            &self,
            _ac: &str,
            start: i32,
            end: Option<i32>,
            _kind: IdentifierType,
        ) -> Result<String, HgvsError> {
            let seq = "ACGTACGTACGTACGTACGT"; // A=0, C=1, G=2, T=3, A=4, ...
            let s = (start.max(0) as usize).min(seq.len());
            let e = end.map_or(seq.len(), |e| (e as usize).min(seq.len()));
            Ok(seq[s..e.max(s)].to_string())
        }
        fn get_symbol_accessions(
            &self,
            _s: &str,
            _f: IdentifierKind,
            _t: IdentifierKind,
        ) -> Result<Vec<(IdentifierType, String)>, HgvsError> {
            Ok(vec![])
        }
        fn get_identifier_type(&self, _id: &str) -> Result<IdentifierType, HgvsError> {
            Ok(IdentifierType::GenomicAccession)
        }
    }

    struct MockSearch;
    impl TranscriptSearch for MockSearch {
        fn get_transcripts_for_region(
            &self,
            _ac: &str,
            _s: i32,
            _e: i32,
        ) -> Result<Vec<String>, HgvsError> {
            Ok(vec![])
        }
    }

    #[test]
    fn normalize_writes_a_repeated_insertion_as_a_duplication() -> Result<(), HgvsError> {
        let hdp = MockDataProvider;
        let mapper = VariantMapper::new(&hdp);
        let norm = |hgvs: &str| -> Result<String, HgvsError> {
            Ok(mapper
                .normalize_variant(crate::parse_hgvs_variant(hgvs)?)?
                .to_string())
        };

        // Reference is ACGT repeated; base 2 (index 1) is C.
        assert_eq!(norm("NC_000001.11:g.2_3insC")?, "NC_000001.11:g.2dupC");
        assert_eq!(norm("NM_000123.4:c.2_3insC")?, "NM_000123.4:c.2dupC");
        assert_eq!(norm("NM_000123.4:n.2_3insC")?, "NM_000123.4:n.2dupC");

        // A whole-unit insertion into a repeat shifts to the 3' end of the run
        // first, then duplicates the last copy.
        assert_eq!(
            norm("NC_000001.11:g.4_5insACGT")?,
            "NC_000001.11:g.17_20dupACGT"
        );

        // An insertion that does not repeat its neighbours stays an insertion.
        assert_eq!(norm("NC_000001.11:g.2_3insTT")?, "NC_000001.11:g.2_3insTT");
        Ok(())
    }

    #[test]
    fn test_normalize_format_question_equals_xaa() {
        let hdp = MockDataProvider;
        let search = MockSearch;
        let eq = VariantEquivalence::new(&hdp, &search);

        // '?' should normalize to 'X', the same as 'Xaa' -> 'X'.
        // This ensures p.Met1? and p.Met1Xaa compare equal after normalization.
        let q = eq.normalize_format("NP_000001.1:p.Met1?");
        let xaa = eq.normalize_format("NP_000001.1:p.Met1Xaa");
        assert_eq!(q, xaa, "p.Met1? and p.Met1Xaa should normalize identically");

        // Parentheses and predicted markers are still stripped.
        let predicted = eq.normalize_format("NP_000001.1:p.(Met1Val)");
        let bare = eq.normalize_format("NP_000001.1:p.Met1Val");
        assert_eq!(predicted, bare);

        // 3-letter codes normalize to 1-letter.
        let three = eq.normalize_format("NP_000001.1:p.Met1Val");
        let one = eq.normalize_format("NP_000001.1:p.M1V");
        assert_eq!(three, one);
    }
}
