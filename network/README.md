##Network

Builds and analyses a directed metabolite interaction network from amino acid
perturbation experiments.

##Structure

```
network/
├── data/
│   └── edgelist_final.txt           # final combined edge list (provided)
│   └── ppr_scores.txt               # PPR influence scores (script output)
├── scripts/
│   ├── NETWORK_CONSTRUCTION.md      # methods record: how edgelist_final.txt was built
│   └── network_influence.R       # PageRank influence scores and plots ← public 
other reference data
└── property-analysis/
    └── data/
        └── AA_properties.csv            # amino acid physicochemical properties
```

How `edgelist_final.txt` was constructed — including all parameters, edge types, and normalisation steps — is documented in `./network/scripts/NETWORK_CONSTRUCTION.md`.
