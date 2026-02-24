## Processing and Statistical Summary

Prepare the coverage BED file needed for analysis.

**Steps:**
1. Upload a **gene list** (.txt, one HGNC symbol per line) or **target BED** (.bed, 4 columns)
2. Upload a **.list** file with one absolute path per line to your BAM or BED coverage files
3. Select reference genome (hg19/hg38) and chromosome notation (chr / number)
4. Select file type — **BAM** (set MAPQ and base quality) or **BED coverage** (set coordinate system)
5. Click **Process Coverage Files**

**Output:** a coverage BED file + statistical summary, both saved in `outDir/output/`. The BED file can be loaded directly into the Coverage Analysis tab.

**After processing:** an optional filter panel lets you subset by sample and threshold and download the result as XLSX.

→ See the full guide for detailed instructions and troubleshooting.
