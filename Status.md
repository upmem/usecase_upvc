# Status

This project is still a work in progress.
This document presents the current work status.

## Current problems

### Mapping order
As read processing is multithreaded, read order is not kept when reporting mappings in sam file.
This makes it somewhat harder to examine sam files as mappings from a same template (read pair) get separated.
However, considering there is a sorting step before the gatk step, this should probably not affect results.

### Mapping score
The mapq we compute is much different from the one bwa computes.
Changing the reported mapq to fit more closely with bwa results would be a good first step.
Looking more in-depth into the reported `AS` value, it could be interesting to change the costs used in dpu computations.
This could change some mappings.

### Unmapped read reporting
Reporting of unmapped reads which are part of a template with another read mapped has been implemented.
(But not properly tested).
Templates with no mappings are not reported at all however.
I am not sure if they should be reported or not.

## Next steps

### GATK evaluation
Once those problems are resolved, generated sam file should be fed back into a gatk pipeline.
Eventual problems which come up need to be resolved.
And ultimately, results (recall & precision) need to be checked.
Then false negatives and false positive calls need to be investigated in order to understand what they are caused by and whether they are caused by fixable issues.
