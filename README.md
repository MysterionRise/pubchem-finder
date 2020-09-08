# pubchem-finder
Tool for finding chemical molecules, compounds, reactions, etc in PubChem (https://pubchem.ncbi.nlm.nih.gov/)

# pubchem-crawler

## Download full pubchem dump

```
python3 pubchem-crawler/crawl.py download --pubchem-dir=/huge/disk
```

## Index pubchem dump into Elasticsearch

```
docker-compose -f elastic/docker-compose.yml up 
python3 pubchem-crawler/crawl.py extract --database=elastic --elastic-no-verify-certs
```