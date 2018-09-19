configfile: 'workflows/config.yaml'

include: 'workflows/shared.snakefile'
include: 'workflows/shared_dm.snakefile'

sys.path.append(os.path.expanduser('~/relmapping/scripts/yarp'))
import yarp as yp

include: 'workflows/tg.snakefile'
include: 'workflows/bwa.snakefile'
include: 'workflows/macs2.snakefile'

include: 'workflows/aln.snakefile'

include: 'workflows/atac.snakefile'
include: 'workflows/atac_HS491.snakefile'
include: 'workflows/atac_cfp1.snakefile'

include: 'workflows/adhoc.snakefile'

#include: 'workflows/dnase_mnase.snakefile'

include: 'workflows/lcap.snakefile'
include: 'workflows/lcap_tm.snakefile'
include: 'workflows/lcap_geo.snakefile'

include: 'workflows/scap.snakefile'
#include: 'workflows/scap_labelled.snakefile'

include: 'workflows/scap_rm_non_coding.snakefile'
include: 'workflows/scap_rm_exonic.snakefile'

include: 'workflows/hmod.snakefile'

#include: 'workflows/seq.snakefile'
#include: 'workflows/dm.snakefile'

include: 'workflows/motifs.snakefile'
include: 'workflows/chip_distribution.snakefile'

include: 'workflows/batch.snakefile'
#include: 'workflows/processed_tracks.snakefile'
include: 'workflows/yapc.snakefile'

include: 'workflows/annot_cb.snakefile'
