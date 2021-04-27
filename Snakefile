import os 
import gzip
import minus80 as m80

#------------------------------------------------------------------------------
#                "Labs grow great when old farts create workflows               
#                 whose ease they know they shall never enjoy."                 
#                                                           -Rob
#------------------------------------------------------------------------------

from src import utils
from pathlib import Path

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id, 
    secret_access_key=s3_access_key,
    keep_local=True
)

configfile: "config.yaml"

ncbi_id_map = {
    'NC_009144.3':'chr1',  'NC_009145.3':'chr2',  'NC_009146.3':'chr3',
    'NC_009147.3':'chr4',  'NC_009148.3':'chr5',  'NC_009149.3':'chr6',
    'NC_009150.3':'chr7',  'NC_009151.3':'chr8',  'NC_009152.3':'chr9',
    'NC_009153.3':'chr10', 'NC_009154.3':'chr11', 'NC_009155.3':'chr12',
    'NC_009156.3':'chr13', 'NC_009157.3':'chr14', 'NC_009158.3':'chr15',
    'NC_009159.3':'chr16', 'NC_009160.3':'chr17', 'NC_009161.3':'chr18',
    'NC_009162.3':'chr19', 'NC_009163.3':'chr20', 'NC_009164.3':'chr21',
    'NC_009165.3':'chr22', 'NC_009166.3':'chr23', 'NC_009167.3':'chr24',
    'NC_009168.3':'chr25', 'NC_009169.3':'chr26', 'NC_009170.3':'chr27',
    'NC_009171.3':'chr28', 'NC_009172.3':'chr29', 'NC_009173.3':'chr30',
    'NC_009174.3':'chr31', 'NC_009175.3':'chrX',  'NC_001640.1':'chrMt'
}

ensem_id_map = utils.gen_ensem_id_map()

# ----------------------------------------------------------
#       It's the all, ya'll
# ----------------------------------------------------------

rule all:
    input:
       # STAR index
        expand('{bucket}/public/refgen/{release}/STAR_INDICES/upload.done',
            bucket=config['bucket'],
            release=[
                'GCF_002863925.1_EquCab3.0',
                'Equus_caballus.EquCab3.0.103'
            ]
        ),
       # Salmon index
        expand('{bucket}/public/refgen/{release}/SALMON_INDEX/upload.done',
            bucket=config['bucket'],
            release=[
                'GCF_002863925.1_EquCab3.0',
                'Equus_caballus.EquCab3.0.103'
            ]
        )

# ----------------------------------------------------------
#       Make "nice" FASTAs
# ----------------------------------------------------------

rule nice_fasta:
    input:
        fna = lambda wildcards: FTP.remote(f'{config[wildcards.release]["fasta"]}')
    output:
        fna_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna.gz')
    threads: 2
    run:
        # get convert dict based on release
        id_map = utils.gen_ensem_id_map()
        if "NCBI" in config[wildcards.release]["resource"]:
            id_map = ncbi_id_map

        with utils.RawFile(input.fna) as IN, \
            gzip.open(output.fna_nice,'wt') as OUT:
            for line in IN:
                if line.startswith('>'):
                    name, *fields = line.lstrip('>').split()
                    if name in id_map:
                        new_name = '>' + id_map[name]
                        line = ' '.join([new_name, name] + fields + ['\n'])
                print(line,file=OUT,end='')

# ----------------------------------------------------------
#       Make "nice" GFFs
# ----------------------------------------------------------

rule nice_gff:
    input:
        gff = lambda wildcards: FTP.remote(f'{config[wildcards.release]["gff"]}')
    output:
        gff_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.gff.gz')
    threads: 2
    run:
        # get convert dict based on release
        id_map = utils.gen_ensem_id_map()
        if 'NCBI' in config[wildcards.release]['resource']:
            id_map = ncbi_id_map
        
        with utils.RawFile(input.gff) as IN, \
            gzip.open(output.gff_nice,'wt') as OUT:
            for line in IN:
                id,*fields = line.split('\t')
                if id in id_map:
                    id = id_map[id]
                print(id,*fields,file=OUT,sep='\t',end='')

# ----------------------------------------------------------
#       Make "nice" transcriptomic FASTAs
# ----------------------------------------------------------

rule make_transcriptomic_fna:
    input:
        fna_nice = S3.remote('{release}/public/refgen/{release}/{release}_genomic.nice.fna.gz'),
        gff_nice = S3.remote('{release}/public/refgen/{release}/{release}_genomic.nice.gff.gz')
    output:
        nice_trx = S3.remote('{bucket}/public/refgen/{release}/{release}_transcriptomic.nice.fna')
    run:
        import locuspocus as lp
        
        # freezable name - get ensembl release number
        freeze = f'ensemblEc3v{wildcards.release.split(".")[-1]}'
        if 'NCBI' in config[wildcards.release]['resource']:
            freeze = 'ncbiEc3'

        # Create the loci db
        fna = lp.Fasta.from_file(freeze,input.fna_nice)
        # Create the GFF db
        gff = lp.Loci(freeze)
        gff.import_gff(input.gff_nice)
        
        # Load the GFF and FNA dbs
        fna = lp.Fasta(freeze) 
        gff = lp.Loci(freeze)
        with open(output.nice_trx,'w') as OUT:
            gff.set_primary_feature_type('gene')
            for gene in gff:
                longest = None
                max_length = 0
                # calulcate the length of each mRNA
                for feature in gene.subloci:
                    # skip non mRNA features
                    if feature.feature_type != 'mRNA':
                        continue
                    # calculate the total length of all the exons that make up the mRNA
                    exon_length = sum([len(x) for x in feature.subloci if x.feature_type == 'exon']) 
                    # Store info it its the longest
                    if exon_length > max_length:
                        longest = feature
                        max_length = exon_length
                if longest is None:
                    continue
                # Print out the nucleotides for the longest mRNA
                print(f">{gene.name}|{feature.name}",file=OUT)
#                exon_seq = ''.join([fna[x.chromosome][x.start:x.end] for x in longest.subloci if x.feature_type == 'exon'])
                exon_seq = ''
                for x in longest.subloci:
                    # add in a print statement to get the bad chromosome
                    print(x.chromosome)
                    if x.feature_type != 'exon':
                        continue
                    exon_seq += fna[x.chromosome][x.start:x.end]
                n = 90
                for chunk in [exon_seq[i:i+n] for i in range(0,len(exon_seq),90)]:
                    print(chunk,file=OUT)
                
# ----------------------------------------------------------
#       Build STAR indices
# ----------------------------------------------------------

rule build_star:
    input:
        fna_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna.gz'),
        gff_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.gff.gz')
    output:
        touch('{bucket}/public/refgen/{release}/STAR_INDICES/upload.done'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/Genome'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/SA'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/SAindex'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/chrLength.txt'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/chrName.txt'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/chrNameLength.txt'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/chrStart.txt'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/exonGeTrInfo.tab'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/exonInfo.tab'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/geneInfo.tab'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/genomeParameters.txt'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/sjdbInfo.txt'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/sjdbList.fromGTF.out.tab'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/sjdbList.out.tab'),
        S3.remote('{bucket}/public/refgen/{release}/STAR_INDICES/transcriptInfo.tab')
    threads: 24
    run:
        # get unzipped fasta and gff paths
        unzip_fna = os.path.splitext(input.fna_nice)[0]
        unzip_gff = os.path.splitext(input.gff_nice)[0]

        shell(f'''
            set -e

            gunzip -c {{input.fna_nice}} > {unzip_fna}
            gunzip -c {{input.gff_nice}} > {unzip_gff}

            STAR \
                --runThreadN {{threads}} \
                --runMode genomeGenerate \
                --genomeDir {{wildcards.bucket}}/public/refgen/{{wildcards.release}}/STAR_INDICES/ \
                --genomeFastaFiles {unzip_fna} \
                --sjdbGTFfile {unzip_gff} \
                --sjdbGTFtagExonParentTranscript Parent
        ''')

# ----------------------------------------------------------
#       Build SALMON indices
# ----------------------------------------------------------

rule prepare_hybrid_fasta:
    input:
        fna_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.fna.gz'),
        gff_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.gff.gz'),
        trx_nice = S3.remote('{bucket}/public/refgen/{release}/{release}_transcriptomic.nice.fna')
    output:
        touch('{bucket}/public/refgen/{release}/SalmonMetadata/upload.done'),
        gtf_nice  = S3.remote('{bucket}/public/refgen/{release}/{release}_genomic.nice.gtf'),
        gentrome  = S3.remote('{bucket}/public/refgen/{release}/SalmonMetadata/gentrome.fa'),
        decoy_ids = S3.remote('{bucket}/public/refgen/{release}/SalmonMetadata/decoys.txt')
    threads: 2
    shell:
        '''
            set -e

            gffread -T {input.gff_nice} \
                -o {output.gtf_nice}

            bash generateDecoyTranscriptome.sh \
                -j {threads} \
                -b /home/.conda/bin/bedtools \
                -m /home/.conda/bin/mashmap \
                -a {output.gtf_nice} \
                -g {input.fna_nice} \
                -t {input.trx_nice} \
                -o {wildcards.bucket}/public/refgen/{wildcards.release}/SalmonMetadata/
        '''

rule build_salmon_index:
    input:
        gentrome  = S3.remote('{bucket}/public/refgen/{release}/SalmonMetadata/gentrome.fa'),
        decoy_ids = S3.remote('{bucket}/public/refgen/{release}/SalmonMetadata/decoys.txt')
    output:
        touch('{bucket}/public/refgen/{release}/SALMON_INDEX/upload.done'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/duplicate_clusters.tsv'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/complete_ref_lens.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/seq.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/rank.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/reflengths.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/ctg_offsets.bin'),
       #S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/eqtable.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/ctable.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/refAccumLengths.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/refseq.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/info.json'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/pos.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/mphf.bin'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/versionInfo.json'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/ref_indexing.log'),
        S3.remote('{bucket}/public/refgen/{release}/SALMON_INDEX/pre_indexing.log')
    threads: 8
    shell:
        '''
            salmon index \
                -i {wildcards.bucket}/public/refgen/{wildcards.release}/SALMON_INDEX \
                -t {input.gentrome} \
                -d {input.decoy_ids} \
                -p {threads}
        '''

