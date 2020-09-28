# pubchem-finder
Tool for finding chemical molecules, compounds, reactions, etc in PubChem (https://pubchem.ncbi.nlm.nih.gov/)

# pubchem-crawler

## Index pubchem dump into Elasticsearch


```
docker-compose -f elastic/docker-compose.yml up 
cd pubchem/
export PYTHONPATH=$(pwd)
python3 pubchem.py pull
```

By default `pubchem` uses `pubchem/` folder in your home directory.
You can change this directory using `--workdir` flag.

`--tmpdir` should be pointed to SSD. If you already use ssd for `--workdir`, 
you can skip this option.   