---
author: Bin He
date: 2024-02-03
---

- YGOB-Pillars.tab is downloaded on 2024-02-02 from http://ygob.ucd.ie/
- CGOB-Pillars.tab is downloaded on 2024-02-02 from http://cgob.ucd.ie/
- CGOB-Pillars-edited.tab is edited (see below)

CGOB-Pillars was edited to correct four lines that have a different number of columns than the rest
```unix
awk '{if (NF != 18) print NR, NF}' CGOB-Pillars.tab
# all other lines have 18 columns, except for the following
17339 0
17340 17
18338 4
18339 13
```

I manually edited the file to merge line 17339 and 17340, adding "---" to the beginning to make that line contain the same number of columns as the others (18 columns). The same was done for 18338 and 18339. The edited file is stored and used for analysis.

Another descrepancy is that `CGOB-README` listed 16 columns (species) while the actual file contained 18. Two, labeled as "PGUG" and "SPAP" were not annotated. This only affects the column selection for us for "S. cerevisiae", which will be col18 instead of col16.
