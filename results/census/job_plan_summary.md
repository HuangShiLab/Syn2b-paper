# GTDB within-species census — task plan

- clusters with >=2 genomes: 22,535
- genomes in those clusters: 274,374
- tasks (query batches of 100): 24,063
- digest budget: 3 core-hours
- synteny budget: 3,556 core-hours
- total: 3,559 core-hours  (at 9 ms/ordered pair)
- capacity per long job: 16 cores x 96 h = 1,536 core-hours
- long jobs needed (including ~20% I/O headroom): 6

## Largest tasks (schedule first; workers claim cost-descending anyway)

| cluster | genomes | task batch | est core-h |
|---|---:|---|---:|
| Escherichia coli | 26,859 | [0,100) | 6.7 |
| Escherichia coli | 26,859 | [100,200) | 6.7 |
| Escherichia coli | 26,859 | [200,300) | 6.7 |
| Escherichia coli | 26,859 | [300,400) | 6.7 |
| Escherichia coli | 26,859 | [400,500) | 6.7 |
| Escherichia coli | 26,859 | [500,600) | 6.7 |
| Escherichia coli | 26,859 | [600,700) | 6.7 |
| Escherichia coli | 26,859 | [700,800) | 6.7 |

Execution policy: submit a FEW LONG array jobs (e.g. `--array=0-5%6`),
each running one `census_worker.py` with `--workers = node cores`. Tasks are
claimed atomically via lock directories on the shared filesystem, so jobs can
start at different times, die and be resubmitted without double work. A space
auditor inside each worker watches quota during the run (see README).
