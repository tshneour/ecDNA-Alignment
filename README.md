## Usage
**In order to use:** Run collect.py to collect reads needed for AA SV refinement. This must be rerun everytime you intend to set the refinement region. The reads will be kept in a file named "alignments.tsv".  
Run refine.py to perform the refinement. The usage of both scripts are shown below. The refinement is currently decided by the candidate end with the most read support. When a tie is encountered, the first is taken.

### collect.py
usage: collect.py [-h] [-v] refine sum bam  

Collect reads from the ends of an AA SV for refinement.  

**Positional Arguments:**
- `refine`: Radius of refinement region
- `sum`: Path to AA amplicon summaries folder
- `bam`: Path to BAM file

**Options:**
- `-h, --help`: Show this help message and exit
- `-v, --verbose`: Include debug output

### refine.py
**Note:** 0 is the default value for verbose. 1 shows read support for the refinement (best) candidate. 2 shows read support for all candidates.  

usage: refine.py [-h] [-l] [-b BREAKPOINTS [BREAKPOINTS ...]] [-v VERBOSE]  

Refine AA SVs.  

**Options:**
- `-h, --help`: Show this help message and exit
- `-l, --list`: Shows all breakpoints
- `-b BREAKPOINTS [BREAKPOINTS ...], --breakpoints BREAKPOINTS [BREAKPOINTS ...]`: List info for these breakpoints (use `--list` to see breakpoint_number)
- `-v VERBOSE, --verbose VERBOSE`: Show information about specific reads