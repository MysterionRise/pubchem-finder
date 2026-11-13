# Batch encoding fixtures

These are synthetic, minimal reproductions of the byte encodings reported in
issue #113, not files exported by Excel.

- `utf8_bom_smiles_first.csv`: UTF-8 BOM, SMILES first, CRLF line endings.
- `utf8_bom_id_first.csv`: UTF-8 BOM, identifier first, CRLF line endings.
- `windows_1252.csv`: no BOM; the name contains a Windows-1252 `0xE9` byte
  (`café`), which is invalid UTF-8.

Keep the original bytes: decoding or re-encoding these files can remove the
condition they reproduce.
