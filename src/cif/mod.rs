mod ast;
mod lexer;
mod parser;

pub use ast::*;

#[cfg(test)]
mod integration_tests {
    use super::*;

    // The block name, `_entry.id` value, and `_atom_site` tag names are
    // checked against the real cifs/4D1M.cif (a pure-protein structure).
    // The `"O5'"` value is a synthetic addition, not something that file
    // contains (it has no primed atoms) — added here to cover the
    // nucleotide-style quoting real PDB entries do exhibit elsewhere.
    // Written as an inline literal, not a disk-read fixture, per preference.
    const ATOM_SITE_SNIPPET: &str = "# minimal PDBx/mmCIF-style block\n\
data_4D1M\n\
_entry.id   4D1M\n\
\n\
loop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.label_atom_id\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
ATOM 1 N      11.104 13.207 2.000\n\
ATOM 2 CA     12.560 13.298 2.500\n\
ATOM 3 \"O5'\" 13.000 14.000 3.000\n\
";

    #[test]
    fn parses_without_panic_and_finds_the_primary_block() {
        let cif = Cif::new(ATOM_SITE_SNIPPET);
        assert!(cif.block("4D1M").is_some());
    }

    #[test]
    fn entry_id_item_is_retrievable() {
        let cif = Cif::new(ATOM_SITE_SNIPPET);
        let block = cif.block("4D1M").unwrap();
        assert!(matches!(block.get("_entry.id"), Some(Value::Str(s)) if s == "4D1M"));
    }

    #[test]
    fn atom_site_column_yields_one_coordinate_per_row() {
        let cif = Cif::new(ATOM_SITE_SNIPPET);
        let block = cif.block("4D1M").unwrap();
        let xs: Vec<f64> = block
            .column("_atom_site.Cartn_x")
            .unwrap()
            .map(|v| v.as_float().unwrap())
            .collect();
        assert_eq!(xs, vec![11.104, 12.560, 13.000]);
    }

    #[test]
    fn minimal_snippet_with_one_item_and_one_loop() {
        let cif = Cif::new("data_m\n_a 1\nloop_\n_b\n_c\n1 2\n3 4");
        let block = cif.block("m").unwrap();
        assert!(matches!(block.get("_a"), Some(Value::Integer { value: 1, .. })));
        let bs: Vec<i64> = block.column("_b").unwrap().map(|v| v.as_int().unwrap()).collect();
        let cs: Vec<i64> = block.column("_c").unwrap().map(|v| v.as_int().unwrap()).collect();
        assert_eq!(bs, vec![1, 3]);
        assert_eq!(cs, vec![2, 4]);
    }

    #[test]
    fn quoted_value_round_trips_through_the_query_api() {
        let cif = Cif::new("data_x\n_a 'has spaces'");
        assert_eq!(cif.block("x").unwrap().get("_a").unwrap().as_str(), Some("has spaces"));
    }

    #[test]
    fn comment_header_and_blank_lines_are_pure_trivia() {
        let plain = Cif::new("data_x\n_a 1");
        let with_trivia = Cif::new("# a header comment\n\n\ndata_x\n_a 1");

        assert_eq!(plain.block("x").unwrap().name, with_trivia.block("x").unwrap().name);
        assert!(matches!(
            (
                plain.block("x").unwrap().get("_a"),
                with_trivia.block("x").unwrap().get("_a")
            ),
            (
                Some(Value::Integer { value: 1, .. }),
                Some(Value::Integer { value: 1, .. })
            )
        ));
    }

    #[test]
    fn dictionary_style_save_frame_still_yields_a_usable_cif() {
        // `SaveFrame` exists in the AST but `try_member` never recognizes one,
        // so parsing silently stops the moment a `save_` heading appears. RCSB
        // *coordinate* mmCIF never contains save frames, but a dictionary file
        // would hit this gap. The member-list shape of the truncation is
        // pinned at the parser unit-test level
        // (parser::tests::save_frame_is_not_parsed_and_silently_ends_the_block);
        // this test instead checks the public-API guarantee: the pipeline
        // doesn't panic and still hands back exactly one usable block, not a
        // corrupted or empty `Cif`.
        let cif = Cif::new("data_x\n_a 1\nsave_s\n_b 2\nsave_\n_c 3");
        assert_eq!(cif.blocks.len(), 1);
        let block = cif.block("x").unwrap();
        assert!(matches!(block.get("_a"), Some(Value::Integer { value: 1, .. })));
    }
}
