# Security policy

Please report security problems privately through
[GitHub's private vulnerability reporting](https://github.com/MysterionRise/pubchem-finder/security/advisories/new),
not as a public issue. You should get a reply within a week.

Fixes are released for the latest version only.

## What PubChem Finder touches

It runs locally. It makes outbound HTTPS requests only to PubChem
(`pubchem.ncbi.nlm.nih.gov`) and, while building an index, to NCBI's download server
(`ftp.ncbi.nlm.nih.gov`). It writes only to its data directory, to files you name with
`--output`, and to the downloads folder inside the data directory.

`pubchem-finder mcp` talks over stdio by default. With `--transport http` it also listens
on a port: `127.0.0.1:8000` unless you choose otherwise. The HTTP server has no
authentication or TLS of its own and warns when bound to anything but loopback; to share
it, put an HTTPS reverse proxy with authentication in front, as
[docs/clients.md](docs/clients.md#over-http) describes.

## How releases are built

Releases are built and published by GitHub Actions from this repository's master branch:
PyPI trusted publishing (no stored tokens) with provenance attestations, a release
environment that needs the maintainer's approval, and actions pinned to commit SHAs.
[docs/releasing.md](docs/releasing.md) has the details.
