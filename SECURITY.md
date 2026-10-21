# Security policy

Please report security problems privately through
[GitHub's private vulnerability reporting](https://github.com/MysterionRise/pubchem-finder/security/advisories/new),
not as a public issue. You should get a reply within a week.

Fixes are released for the latest version only.

PubChem Finder runs locally. It makes outbound HTTPS requests only to PubChem
(`pubchem.ncbi.nlm.nih.gov`) and, while building an index, to NCBI's download server
(`ftp.ncbi.nlm.nih.gov`). It writes only to its data directory, to files you name with
`--output`, and to the downloads folder inside the data directory.
