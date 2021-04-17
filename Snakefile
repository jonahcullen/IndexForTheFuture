container: "docker://jonahcullen/ec3index:v0.1.2"

localrules: 
    nice_ncbi_fasta,nice_ensembl_fasta, \
    nice_ncbi_gff,nice_ensembl_gff, \
    make_transcriptomic_fna_ncbi,make_locpoc_dbs_ncbi, \
    make_transcriptomic_fna_ensembl,make_locpoc_dbs_ensembl

import os 
#import gzip
#import minus80 as m80

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
    secret_access_key=s3_access_key
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

# refseq - GCF v genbank - GCA
#NCBI_ASSEM = config['NCBI']['ASSEMBLY']
#ENSEMBL_ASSEM = config['ENSEMBL']['ASSEMBLY']
#BUCKET = config['BUCKET']
#
#ASSEM=[config['NCBI']['ASSEMBLY'],config['ENSEMBL']['ASSEMBLY']]
#print(ASSEM)

#print(config["Ensembl"]["fasta"])
#print(config['NCBI']['fasta'])



# ----------------------------------------------------------
#       It's the all, ya'll
# ----------------------------------------------------------

rule all:
    input:
       #S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_transcriptomic.nice.fna'),
       #S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_transcriptomic.nice.fna'),
       #S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/gentrome.fa')
       #f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/download.done'
       ## PREVIOUS WORKING ALL INPUT
       #f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/download.done',
       #f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/download.done'
        # NICE FASTA
       #S3.remote(
       #    expand(
       #        '{bucket}/public/refgen/{release}/{release}_genomic.nice.fna.gz',
       #        bucket=config['bucket'],
       #        #resource=["Ensembl","NCBI"],
       #        release=[
       #            "GCF_002863925.1_EquCab3.0",
       #            "Equus_caballus.EquCab3.0.103"
       #        ]
       #    ),
       #    keep_local=True
       #),
       ## NICE GFF
       #S3.remote(
       #    expand(
       #        '{bucket}/public/refgen/{release}/{release}_genomic.nice.gff.gz',
       #        bucket=config['bucket'],
       #        #resource=["Ensembl","NCBI"],
       #        release=[
       #            "GCF_002863925.1_EquCab3.0",
       #            "Equus_caballus.EquCab3.0.103"
       #        ]
       #    ),
       #    keep_local=True
       #),
       ## LOCUSPOCUS DBS
        expand(
            '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Fasta.{release}/thawed/tinydb.json',
            bucket=config['bucket'],
            release=[
                "GCF_002863925.1_EquCab3.0",
                "Equus_caballus.EquCab3.0.103"
            ]
        )
       ## NICE TRANSCRIPTOMIC
       #S3.remote(
       #    expand(
       #        '{bucket}/public/refgen/{release}/{release}_transcriptomic.nice.fna',
       #        bucket=config['bucket'],
       #        release=[
       #            "GCF_002863925.1_EquCab3.0",
       #            "Equus_caballus.EquCab3.0.103"
       #        ]
       #    ),
       #    keep_local=True
       #)

       #expand('{bucket}/public/refgen/{u.resource}/{u.release}/{u.release}_genomic.nice.fna.gz',
       #    u=units.itertuples(),
       #    bucket=config['bucket']
       #)

# ----------------------------------------------------------
#       Make "nice" FASTAs
# ----------------------------------------------------------

rule nice_fasta:
    input:
       #fna = FTP.remote(config['NCBI']['fasta']),
       #dumb = FTP.remote(config['Ensembl']['fasta'])
       #fna = FTP.remote(*get_fastq())
       #fna = FTP.remote(f'{{config["ASSEM"]["FASTA"]}}')
       #fna = FTP.remote(config["ASSEM"]["FASTA"])
       #fna = FTP.remote(f'{{config[{wildcards.release}]["fasta"]}}')
        fna = lambda wildcards: FTP.remote(f'{config[wildcards.release]["fasta"]}')
       #lambda wildcards: \
       #    [f"{config[res]['fasta']}" for res in resource]
       #fna = FTP.remote(f'{config[resource]["fasta"]}')
       #fna = lambda wildcards: FTP.remote("{SLUG}", SLUG=config[wildcards.resource]['fasta'])
       #fna = lambda wildcards: FTP.remote(config[rel]['fasta'] for rel in rsrc_rel)
    output:
       #fna = S3.remote(f'{BUCKET}/public/refgen/{{ASSEM}}/{{ASSEM}}_genomic.nice.fna.gz',keep_local=True)
       #fna = S3.remote('{bucket}/public/refgen/{resource}/{release}/{release}_genomic.nice.fna.gz',keep_local=True)
        fna_nice = S3.remote(
            '{bucket}/public/refgen/{release}/{release}_genomic.nice.fna.gz'
        )
    threads: 2
    resources:
        time   = 120,
        mem_mb = 6000
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

# cram = lambda wildcards: FTP.remote(expand("ftp://ftp.sra.ebi.ac.uk/vol1/{cram}", cram = config['Samples'][wildcards.sample]['WES'])),

#rule get_patient_output:
#        input:
#            lambda wildcards: \
        #                ["mutation_data/{0}_output.tsv".format(gene_id) \
        #                    for gene_id in patient_gene_mutations[wildcards.patient_id]
#                ]
#        output:
#            "patient_data/{patient_id}_output.tsv"
#        shell:
#            """
#            touch {output}
#                    """


# ----------------------------------------------------------
#       Make "nice" GFFs
# ----------------------------------------------------------

rule nice_gff:
    input:
       #gff = FTP.remote(f'{config["NCBI"]["GFF"]}')
        gff = lambda wildcards: FTP.remote(f'{config[wildcards.release]["gff"]}')
    output:
       #gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz',keep_local=True)
        gff_nice = S3.remote(
            '{bucket}/public/refgen/{release}/{release}_genomic.nice.gff.gz'
        )
    threads: 2
    resources:
        time   = 120,
        mem_mb = 6000
    run:
        # get convert dict based on release
        id_map = utils.gen_ensem_id_map()
        if "NCBI" in config[wildcards.release]["resource"]:
            id_map = ncbi_id_map
        
        with utils.RawFile(input.gff) as IN, \
            gzip.open(output.gff_nice,'wt') as OUT:
            for line in IN:
                id,*fields = line.split('\t')
                if id in id_map:
                    id = id_map[id]
                print(id,*fields,file=OUT,sep='\t',end='')


#rule nice_ensembl_gff:
#    input:
#        gff = FTP.remote(f'{config["ENSEMBL"]["GFF"]}')
#    output:
#        gff = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff.gz',keep_local=True)
#    run:
#        with utils.RawFile(input.gff) as IN, \
#            gzip.open(output.gff,'wt') as OUT:
#            for line in IN:
#                id,*fields = line.split('\t')
#                if id in ensem_id_map:
#                    id = ensem_id_map[id]
#                print(id,*fields,file=OUT,sep='\t',end='')
#
## ----------------------------------------------------------
##       Make the LocusPocus Databases for the GFF/Fasta (NCBI)
## ----------------------------------------------------------

rule make_locpoc_dbs_ncbi:
    input:
       #fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz',keep_local=True),
       #gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz',keep_local=True)
        fna_nice = S3.remote(
            f'{config["bucket"]}/public/refgen/{{release}}/{{release}}_genomic.nice.fna.gz'
        ),
        gff_nice = S3.remote(
            f'{config["bucket"]}/public/refgen/{{release}}/{{release}}_genomic.nice.gff.gz'
        )
    output:
       #fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ncbiEquCab3/thawed/tinydb.json',
       #gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ncbiEquCab3/thawed/tinydb.json'
        fna = '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Loci.{release}/thawed/tinydb.json',
        gff = '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Fasta.{release}/thawed/tinydb.json'
       #fna = '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Loci.{resource}/thawed/tinydb.json',
       #gff = '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Fasta.{resource}/thawed/tinydb.json'
    run:
        # Create the loci db
        import locuspocus as lp
        fna = lp.Fasta.from_file('ncbiEquCab3', input.fna)  
        # Create the GFF db
        gff = lp.Loci('ncbiEquCab3')
        gff.import_gff(input.gff)

rule make_transcriptomic_fna_ncbi:
    input:
        #fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ncbiEquCab3/thawed/tinydb.json',
        #gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ncbiEquCab3/thawed/tinydb.json'
         fna = '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Loci.{release}/thawed/tinydb.json',
         gff = '/panfs/roc/groups/0/fried255/cull0084/.minus80/datasets/v1/Fasta.{release}/thawed/tinydb.json'
    output:
         fna = S3.remote('{bucket}/public/refgen/{release}/{release}_transcriptomic.nice.fna',keep_local=True)
#    run:
#        import locuspocus as lp
#        # Load the GFF and FNA dbs
#        fna = lp.Fasta('ncbiEquCab3') 
#        gff = lp.Loci('ncbiEquCab3')
#        with open(output.fna,'w') as OUT:
#            gff.set_primary_feature_type('gene')
#            for gene in gff:
#                longest = None
#                max_length = 0
#                # calulcate the length of each mRNA
#                for feature in gene.subloci:
#                    # skip non mRNA features
#                    if feature.feature_type != 'mRNA':
#                        continue
#                    # calculate the total length of all the exons that make up the mRNA
#                    exon_length = sum([len(x) for x in feature.subloci if x.feature_type == 'exon']) 
#                    # Store info it its the longest
#                    if exon_length > max_length:
#                        longest = feature
#                        max_length = exon_length
#                if longest is None:
#                    continue
#                # Print out the nucleotides for the longest mRNA
#                print(f">{gene.name}|{feature.name}",file=OUT)
##                exon_seq = ''.join([fna[x.chromosome][x.start:x.end] for x in longest.subloci if x.feature_type == 'exon'])
#                exon_seq = ''
#                for x in longest.subloci:
#                    # add in a print statement to get the bad chromosome
#                    print(x.chromosome)
#                    if x.feature_type != 'exon':
#                        continue
#                    exon_seq += fna[x.chromosome][x.start:x.end]
#                n = 90
#                for chunk in [exon_seq[i:i+n] for i in range(0,len(exon_seq),90)]:
#                    print(chunk,file=OUT)
#                
#rule make_locpoc_dbs_ncbi:
#    input:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz',keep_local=True),
#        gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz',keep_local=True)
#    output:
#        fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ncbiEquCab3/thawed/tinydb.json',
#        gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ncbiEquCab3/thawed/tinydb.json'
#    run:
#        # Create the loci db
#        import locuspocus as lp
#        fna = lp.Fasta.from_file('ncbiEquCab3', input.fna)  
#        # Create the GFF db
#        gff = lp.Loci('ncbiEquCab3')
#        gff.import_gff(input.gff)
#        
## ----------------------------------------------------------
##       Make the LocusPocus Databases for the GFF/Fasta (ENSEMBL)
## ----------------------------------------------------------
#
#rule make_transcriptomic_fna_ensembl:
#    input:
#        fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ensemblEquCab3/thawed/tinydb.json',
#        gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ensemblEquCab3/thawed/tinydb.json'
#    output:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_transcriptomic.nice.fna',keep_local=True)
#    run:
#        import locuspocus as lp
#        # Load the GFF and FNA dbs
#        fna = lp.Fasta('ensemblEquCab3') 
#        gff = lp.Loci('ensemblEquCab3')
#        with open(output.fna,'w') as OUT:
#            gff.set_primary_feature_type('gene')
#            for gene in gff:
#                longest = None
#                max_length = 0
#                # calulcate the length of each mRNA
#                for feature in gene.subloci:
#                    # skip non mRNA features
#                    if feature.feature_type != 'mRNA':
#                        continue
#                    # calculate the total length of all the exons that make up the mRNA
#                    exon_length = sum([len(x) for x in feature.subloci if x.feature_type == 'exon']) 
#                    # Store info it its the longest
#                    if exon_length > max_length:
#                        longest = feature
#                        max_length = exon_length
#                if longest is None:
#                    continue
#                # Print out the nucleotides for the longest mRNA
#                print(f">{gene.name}|{feature.name}",file=OUT)
##                exon_seq = ''.join([fna[x.chromosome][x.start:x.end] for x in longest.subloci if x.feature_type == 'exon'])
#                exon_seq = ''
#                for x in longest.subloci:
#                    # add in a print statement to get the bad chromosome
#                    print(x.chromosome)
#                    if x.feature_type != 'exon':
#                        continue
#                    exon_seq += fna[x.chromosome][x.start:x.end]
#                n = 90
#                for chunk in [exon_seq[i:i+n] for i in range(0,len(exon_seq),90)]:
#                    print(chunk,file=OUT)
#
#rule make_locpoc_dbs_ensembl:
#    input:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna.gz',keep_local=True),
#        gff = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff.gz',keep_local=True)
#    output:
#        fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ensemblEquCab3/thawed/tinydb.json',
#        gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ensemblEquCab3/thawed/tinydb.json'
#    run:
#        # Create the loci db
#        import locuspocus as lp
#        fna = lp.Fasta.from_file('ensemblEquCab3', input.fna)  
#        # Create the GFF db
#        gff = lp.Loci('ensemblEquCab3')
#        gff.import_gff(input.gff)
#
## ----------------------------------------------------------
##       Build STAR indices
## ----------------------------------------------------------
#
#rule build_star_ncbi:
#    input:
#        gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz'),
#        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz')
#    output:
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/Genome'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/SA'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/SAindex'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrLength.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrName.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrNameLength.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrStart.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/exonGeTrInfo.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/exonInfo.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/geneInfo.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/genomeParameters.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/sjdbInfo.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/sjdbList.fromGTF.out.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/sjdbList.out.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/transcriptInfo.tab')
#    threads: int(f'{int(config["THREADS"]["STAR"])}')
#    shell:
#        f'''
#          STAR \
#            --runThreadN {{threads}} \
#            --runMode genomeGenerate \
#            --genomeDir {BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/ \
#            --genomeFastaFiles {{input.fna}} \
#            --sjdbGTFfile {{input.gff}} \
#            --sjdbGTFtagExonParentTranscript Parent
#        '''
#
#rule build_star_ensembl:
#    input:
#        gff = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff3.gz'),
#        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna.gz')
#    output:
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/Genome'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/SA'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/SAindex'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrLength.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrName.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrNameLength.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrStart.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/exonGeTrInfo.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/exonInfo.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/geneInfo.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/genomeParameters.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/sjdbInfo.txt'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/sjdbList.fromGTF.out.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/sjdbList.out.tab'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/transcriptInfo.tab') 
#    threads: int(f'{int(config["THREADS"]["STAR"])}')
#    shell:
#        f'''
#          STAR \
#            --runThreadN {{threads}} \
#            --runMode genomeGenerate \
#            --genomeDir {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/ \
#            --genomeFastaFiles {{input.fna}} \
#            --sjdbGTFfile {{input.gff}} \
#            --sjdbGTFtagExonParentTranscript Parent 
#        '''
#
## ----------------------------------------------------------
##       Build SALMON indices (NCBI)
## ----------------------------------------------------------
#
#rule prepare_hybrid_fasta_ncbi:
#    input:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz'),
#        trx = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_transcriptomic.nice.fna'),
#        gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz')
#    output:
#        touch(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/download.done'),
#        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/gentrome.fa'),
#        decoy_ids = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/decoys.txt')
#    threads: int(f'{int(config["THREADS"]["SALMON"])}')
#    shell:
#        f'''
#          gunzip -c {{input.fna}} > {BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna &&
#          gunzip -c {{input.gff}} > {BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff &&
#          gffread -T {BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff -o {BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gtf &&
#          bash generateDecoyTranscriptome.sh \
#            -j {{threads}} \
#            -b /home/.conda/bin/bedtools \
#            -m /home/.conda/bin/mashmap \
#            -a {BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gtf \
#            -g {BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna \
#            -t {{input.trx}} \
#            -o {BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/
#        '''
#
#rule build_salmon_ncbi:
#    input:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/gentrome.fa'),
#        decoy_ids = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SalmonMetadata/decoys.txt')
#    output:
#        touch(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/download.done'),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/duplicate_clusters.tsv',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/complete_ref_lens.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/seq.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/rank.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/reflengths.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/ctg_offsets.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/eqtable.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/ctable.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/refAccumLengths.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/refseq.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/info.json',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/pos.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/mphf.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/versionInfo.json',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/ref_indexing.log',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX/pre_indexing.log',keep_local=True)
#    threads: int(f'{int(config["THREADS"]["SALMON"])}')
#    shell:
#        f'''
#          salmon index \
#            -i {BUCKET}/public/refgen/{NCBI_ASSEM}/SALMON_INDEX \
#            -t {{input.fna}} \
#            -d {{input.decoy_ids}} \
#            -p {{threads}}
#         '''
#
## ----------------------------------------------------------
##       Build SALMON indices (ENSEMBL)
## ----------------------------------------------------------
#
#rule prepare_hybrid_fasta_ensembl:
#    input:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna.gz'),
#        trx = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_transcriptomic.nice.fna'),
#        gff = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff.gz')
#    output:
#        touch(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SalmonMetadata/download.done'),
#        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SalmonMetadata/gentrome.fa'),
#        decoy_ids = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SalmonMetadata/decoys.txt')
#    threads: int(f'{int(config["THREADS"]["SALMON"])}')
#    shell:
#        f'''
#          gunzip -c {{input.fna}} > {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna &&
#          gunzip -c {{input.gff}} > {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff &&
#          gffread -T {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff -o {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gtf &&
#          bash generateDecoyTranscriptome.sh \
#            -j {{threads}} \
#            -b /home/.conda/bin/bedtools \
#            -m /home/.conda/bin/mashmap \
#            -a {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gtf \
#            -g {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna \
#            -t {{input.trx}} \
#            -o {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SalmonMetadata/
#        '''
#
#rule build_salmon_ensembl:
#    input:
#        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SalmonMetadata/gentrome.fa'),
#        decoy_ids = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SalmonMetadata/decoys.txt')
#    output:
#        touch(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/download.done'),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/duplicate_clusters.tsv',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/complete_ref_lens.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/seq.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/rank.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/reflengths.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/ctg_offsets.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/eqtable.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/ctable.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/refAccumLengths.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/refseq.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/info.json',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/pos.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/mphf.bin',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/versionInfo.json',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/ref_indexing.log',keep_local=True),
#        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX/pre_indexing.log',keep_local=True)
#    threads: int(f'{int(config["THREADS"]["SALMON"])}')
#    shell:
#        f'''
#          salmon index \
#            -i {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/SALMON_INDEX \
#            -t {{input.fna}} \
#            -d {{input.decoy_ids}} \
#            -p {{threads}}
#         '''
#
